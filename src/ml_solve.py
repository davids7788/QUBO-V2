import numpy as np
import matplotlib.pyplot as plt
import time

from pattern.doublet import Doublet
from pattern.triplet import Triplet

from tensorflow.keras import optimizers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import model_from_json
from tensorflow.keras import backend

from sklearn.metrics import roc_curve, auc


# nicely formatted time string
def hms_string(sec_elapsed):
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)


train_0 = np.load("../ml_folders/triplet_list_train_0.npy", allow_pickle=True)
train_1 = np.load("../ml_folders/triplet_list_train_0.npy", allow_pickle=True)
train_2 = np.load("../ml_folders/triplet_list_train_0.npy", allow_pickle=True)

t_samples = [train_0, train_1, train_2]

val = np.load("../ml_folders/triplet_list_val.npy", allow_pickle=True)

test = np.load("../ml_folders/triplet_list_test.npy", allow_pickle=True)

training_sample = []
training_truth = []

val_sample = []
val_truth = []

test_sample = []
test_truth = []

test_CWoLa_feature = []
test_truth_CWoLa_feature = []


for sample in t_samples:
    for triplet in sample:
        interactions = list(triplet.interactions.values())
        if not interactions:
            interactions.append(0)
        if not [x for x in interactions if x < 0.0]:
            interactions.append(-1e10)

        feature_vector = [triplet.doublet_1.hit_1_position[0],
                          min([x for x in interactions if x < 0]),
                          max([x for x in interactions if x < 0]),
                          len([x for x in interactions if x < 0]),
                          list(triplet.interactions.values()).count(1)]
        training_sample.append(feature_vector)
        if triplet.quality > 0:
            training_truth.append(0)
        else:
            training_truth.append(1)

for triplet in val:
    interactions = list(triplet.interactions.values())
    if not interactions:
        interactions.append(0)
    if not [x for x in interactions if x < 0.0]:
        interactions.append(-1e10)
    feature_vector = [triplet.doublet_1.hit_1_position[0],
                      min([x for x in interactions if x < 0]),
                      max([x for x in interactions if x < 0]),
                      len([x for x in interactions if x < 0]),
                      list(triplet.interactions.values()).count(1)]
    val_sample.append(feature_vector)
    if triplet.is_correct_match:
        val_truth.append(1)
    else:
        val_truth.append(0)

for triplet in test:
    interactions = list(triplet.interactions.values())
    if not interactions:
        interactions.append(0)
    if not [x for x in interactions if x < 0.0]:
        interactions.append(-1e10)
    feature_vector = [triplet.doublet_1.hit_1_position[0],
                      min([x for x in interactions if x < 0]),
                      max([x for x in interactions if x < 0]),
                      len([x for x in interactions if x < 0]),
                      list(triplet.interactions.values()).count(1)]
    test_sample.append(feature_vector)
    if triplet.is_correct_match:
        test_truth.append(1)
    else:
        test_truth.append(0)
    if triplet.quality > 0:
        test_truth_CWoLa_feature.append(0)
    else:
        test_truth_CWoLa_feature.append(1)


def create_model(optimizer,
                 nodes_list=(16, 32, 64, 32, 16, 8),
                 batch_normalization=True,
                 activation_function='selu',
                 kernel_initializer='glorot_normal',
                 bias_initializer='glorot_normal',
                 batch_size=1024,
                 loss='sparse_categorical_crossentropy',
                 metrics='accuracy'):

    """ creates a sequential model and writes a file with some of the training history
    :param nodes_list list of nodes for each dense layer, has to be the size of num_dense_layers
    :param batch_normalization if True batch normalization is used after each layer
    :param activation_function the used activation function for all but the last layer
    :param kernel_initializer kernel initializer, default 'he_uniform' like the one in the CWoLa paper
    :param bias_initializer bias initializer, default 'he_uniform' like the one in the CWoLa paper
    :param batch_size size of the used batches for training
    :param optimizer self build optimizer or adam with default values
    :param loss function of the model
    :param metrics tuple of metrics
    :return created model"""

    model = Sequential()
    for i in range(len(nodes_list)):
        model.add(Dense(nodes_list[i],
                        activation=activation_function,
                        kernel_initializer=kernel_initializer,
                        bias_initializer=bias_initializer,
                        batch_size=batch_size))
        if batch_normalization and i != len(nodes_list) - 1:      # Never add batch normalization before softmax layer
            model.add(BatchNormalization())

    model.add(Dense(2, activation='softmax'))
    model.compile(optimizer=optimizer, loss=loss, metrics=[metrics])

    return model


# specifying parameters for optimizer
LR = 1e-3
EPOCHS = 75
DECAY = LR /EPOCHS     # in some cases LR / EPOCHS might be useful
BETA_1 = 0.8
BETA_2 = 0.99
AMSGRAD = True

# Preparing k-fold cross validation
K_FOLDING = 5 # optimal value might be 10, but this needs a lot of time

adam = optimizers.Adam(lr=LR, decay=DECAY, beta_1=BETA_1, beta_2=BETA_2, amsgrad=AMSGRAD)

print('############################################################')
print('Using adam optimizer with the following parameters:\n' +
      'learning rate:', LR, '\n'
      'decay:', DECAY,  '\n'
      'beta1:', BETA_1, '\n'
      'beta2:', BETA_2, '\n'
      'amsgrad:', AMSGRAD)
print('############################################################', '\n')


# model loss

file_path_val_loss = 'weights_best_val_loss.h5'
file_path_val_accuracy = 'weights_best_val_accuracy.h5'

start = time.time()
classifier = create_model(adam,
                          nodes_list=(16, 32, 64, 32, 16, 8),
                          batch_normalization=True,
                          activation_function='selu',
                          kernel_initializer='glorot_normal',
                          bias_initializer='glorot_normal',
                          batch_size=256,
                          loss='sparse_categorical_crossentropy',
                          metrics='accuracy')

# checking for lowest val loss and accuracy
checkpoint_val_loss = \
    ModelCheckpoint(file_path_val_loss,
                    monitor='val_loss',
                    verbose=1,
                    save_best_only=True,
                    mode='min')

checkpoint_val_accuracy = \
    ModelCheckpoint(file_path_val_accuracy,
                    monitor='val_accuracy',
                    verbose=1,
                    save_best_only=True,
                    mode='max')

callbacks_list = [checkpoint_val_loss, checkpoint_val_accuracy]
model_history = classifier.fit(training_sample,
                               training_truth,
                               class_weight={0: training_truth.count(1),
                                             1: training_truth.count(0)},
                               epochs=EPOCHS,
                               verbose=2,
                               validation_data=(val_sample,
                                                val_truth),
                               callbacks=callbacks_list,
                               shuffle=True)

end = time.time()

model_json = classifier.to_json()
with open('model.json', 'w') as json_file:
    json_file.write(model_json)
classifier.save_weights('after_last_epoch.h5')


# Plot training & validation accuracy values
plt.figure()
plt.plot(model_history.history['accuracy'])
plt.plot(model_history.history['val_accuracy'])
plt.title('Model accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='best')
plt.savefig('model_history-model accuracy')
plt.show()

# Plot training & validation loss values
plt.figure()
plt.plot(model_history.history['loss'])
plt.plot(model_history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='best')
plt.savefig('model_history-model accuracy')
plt.show()

print(f'elapsed time: {hms_string(end - start)}')

json_file = open('model.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)

# load weights into new model

loaded_model.load_weights('weights_best_val_loss.h5')
predictions = loaded_model.predict(test_sample)
fpr, tpr, thresholds = roc_curve(test_truth,
                                 predictions[:, 1])
roc_auc = auc(fpr, tpr)

predictions_CWoLa_feature = loaded_model.predict(test_sample)
fpr_CWoLa_feature, tpr_CWoLa_feature, thresholds_CWoLa_feature = \
    roc_curve(test_truth_CWoLa_feature,  predictions_CWoLa_feature[:, 1])
roc_auc_CWoLa_feature = auc(fpr_CWoLa_feature, tpr_CWoLa_feature)

plt.figure()
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
plt.xlim([0, 1.0])
plt.ylim([0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC-Curve')
plt.legend(loc="lower right")
plt.savefig('ROC_Curve')


# Training data

plt.plot(fpr_CWoLa_feature, tpr_CWoLa_feature, color='b',
         label=r'CWoLa feature separation (AUC = %0.2f)' % roc_auc_CWoLa_feature,
         lw=2, alpha=.8)

# Evaluation data

plt.plot(fpr, tpr, color='g',
         label=r'Test S/BG (AUC = %0.2f)' % roc_auc,
         lw=2, alpha=.8)

plt.legend(loc="best")
plt.show()
