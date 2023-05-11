import numpy as np
import matplotlib.pyplot as plt
import time
import os


from pattern.doublet import Doublet
from pattern.triplet import Triplet
from utility.time_tracking import hms_string

from tensorflow.keras import optimizers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import model_from_json
from tensorflow.keras import backend
from track_reconstruction.create_reco_xplets import reco_xplets_simplified_LUXE
from track_reconstruction.track_reconstruction_efficiency import track_reconstruction_efficiency_simplified_LUXE

from sklearn.metrics import roc_curve, auc

# train_data_list = ["10", "11", "12", "13", "14", "15", "16", "17"]
# val_data_list = ["18", "19"]
# test_data_list = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09"]

xi = 5.0
train_data_list = ["10", "11", "12"]
val_data_list = ["18"]
test_data_list = ["00"]

directory = "../../qubo_ml_datasets/xi_5"
files_in_directory = os.listdir(directory)

new_folder = f"{directory}/ml_solve_{str(np.random.randint(1e8, 1e9))}"
print(f"Creating folder: {new_folder}\n")
os.mkdir(new_folder)

print(f"Starting ml_solve for xi={xi}:")

training_samples = []
validation_samples = []
# test_samples = ["../../qubo_ml_datasets/xi_7/e0gpc_7.0_0000_sl-example/triplet_list.npy"]
test_samples = []


for folder in files_in_directory:
    if not os.path.isdir(f"{directory}/{folder}"):
        continue
    if folder.split("_")[2][2:] in train_data_list:
        training_samples.append(f"{directory}/{folder}/triplet_list.npy")
    if folder.split("_")[2][2:] in val_data_list:
        validation_samples.append(f"{directory}/{folder}/triplet_list.npy")
    if folder.split("_")[2][2:] in test_data_list:
        test_samples.append(f"{directory}/{folder}/triplet_list.npy")

print("\nSample(s) used for training:")
for sample in training_samples:
    print(sample)
print("\nSample(s) used for validation:")
for sample in validation_samples:
    print(sample)
print("\nSample(s) used for validation:")
for sample in test_samples:
    print(sample)

# norm training samples
conflicts_train = []
connections_train = []
for train_sample in training_samples:
    train_file = np.load(train_sample, allow_pickle=True)
    for triplet in train_file:
        interactions = list(triplet.interactions.values())
        connections = 0
        conflicts = 0
        for value in interactions:
            if value < 0:
                connections += 1
            else:
                conflicts += 1
        conflicts_train.append(conflicts)
        connections_train.append(connections)

max_num_conflicts_train = max(conflicts_train)
max_num_connections_train = max(connections_train)

# norm test samples
conflicts_test = []
connections_test = []
for test_sample in test_samples:
    test_file = np.load(test_sample, allow_pickle=True)
    for triplet in test_file:
        interactions = list(triplet.interactions.values())
        connections = 0
        conflicts = 0
        for value in interactions:
            if value < 0:
                connections += 1
            else:
                conflicts += 1
        conflicts_test.append(conflicts)
        connections_test.append(connections)

max_num_conflicts_test = max(conflicts_test)
max_num_connections_test = max(connections_test)


def extract_feature_vector(triplet, max_conflicts, max_connections):
    interactions = list(triplet.interactions.values())
    connections = [0]
    conflicts = [1]
    for v in interactions:
        if v < 0:
            connections.append(v)
        else:
            conflicts.append(v)
    quality = triplet.quality

    return [triplet.doublet_1.hit_1_position[0],
            triplet.doublet_1.hit_1_position[1],
            100 * quality,
            min(connections),
            max(connections),
            len(connections) / max_connections,
            len(conflicts) / max_conflicts]


train_input_labels = []
train_input_feature_vector = []

val_input_labels = []
val_input_feature_vector = []

test_input_labels = []
test_input_feature_vector = []

test_triplets_complete = []


for train_sample in training_samples:
    train_file = np.load(train_sample, allow_pickle=True)
    for t in train_file:
        train_input_feature_vector.append(extract_feature_vector(t,
                                                                 max_num_conflicts_train,
                                                                 max_num_connections_train))
        if t.is_correct_match():
            train_input_labels.append(1)
        else:
            train_input_labels.append(0)

for val_sample in validation_samples:
    val_file = np.load(val_sample, allow_pickle=True)
    for t in val_file:
        val_input_feature_vector.append(extract_feature_vector(t,
                                                               max_num_conflicts_train,
                                                               max_num_connections_train))
        if t.is_correct_match():
            val_input_labels.append(1)
        else:
            val_input_labels.append(0)

for test_sample in test_samples:
    test_file = np.load(test_sample, allow_pickle=True)
    for t in test_file:
        test_input_feature_vector.append(extract_feature_vector(t,
                                                                max_num_conflicts_test,
                                                                max_num_connections_test))
        test_triplets_complete.append(t)
        if t.is_correct_match():
            test_input_labels.append(1)
        else:
            test_input_labels.append(0)


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
EPOCHS = 100
DECAY = LR / EPOCHS     # in some cases LR / EPOCHS might be useful
BETA_1 = 0.8
BETA_2 = 0.99
AMSGRAD = True


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

file_path_val_loss = f"{new_folder}/weights_best_val_loss.h5"
file_path_val_accuracy = f"{new_folder}/weights_best_val_accuracy.h5"

start = time.time()
classifier = create_model(adam,
                          nodes_list=(8, 16, 32, 32, 16, 8),
                          batch_normalization=True,
                          activation_function='selu',
                          kernel_initializer='glorot_normal',
                          bias_initializer='glorot_normal',
                          batch_size=1024,
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
model_history = classifier.fit(np.asarray(train_input_feature_vector).astype("float32"),
                               np.asarray(train_input_labels).astype("float32"),
                               class_weight={0: train_input_labels.count(1) / len(train_input_labels),
                                             1: train_input_labels.count(0) / len(train_input_labels)},
                               epochs=EPOCHS,
                               verbose=2,
                               validation_data=(np.asarray(val_input_feature_vector).astype("float32"),
                                                np.asarray(val_input_labels).astype("float32")),
                               callbacks=callbacks_list,
                               shuffle=True)

end = time.time()

model_json = classifier.to_json()
with open(f"{new_folder}/model.json", 'w') as json_file:
    json_file.write(model_json)
classifier.save_weights(f"{new_folder}/after_last_epoch.h5")


# Plot training & validation accuracy values
plt.figure()
plt.plot(model_history.history['accuracy'])
plt.plot(model_history.history['val_accuracy'])
plt.title('Model accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='best')
plt.savefig(f"{new_folder}/model_history-model accuracy")
# plt.show()

# Plot training & validation loss values
plt.figure()
plt.plot(model_history.history['loss'])
plt.plot(model_history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Val'], loc='best')
plt.savefig(f"{new_folder}/model_history-model accuracy")
# plt.show()

print(f'elapsed time: {hms_string(end - start)}')

json_file = open(f"{new_folder}/model.json", "r")
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)

# load weights into new model

loaded_model.load_weights(f"{new_folder}/weights_best_val_loss.h5")
predictions = loaded_model.predict(test_input_feature_vector)
fpr, tpr, thresholds = roc_curve(test_input_labels,
                                 predictions[:, 1])
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
plt.xlim([0, 1.0])
plt.ylim([0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC-Curve')
# Evaluation data

plt.plot(fpr, tpr, color='g',
         label=r'Test S/BG (AUC = %0.2f)' % roc_auc,
         lw=2, alpha=.8)

plt.legend(loc="best")
plt.savefig(f"{new_folder}/ROC_Curve")
# plt.show()

# predictions
plt.figure()
plt.hist(predictions[:, 1], bins=50)
plt.xlabel('Prediction')
plt.ylabel('Counts')
plt.legend(loc="lower right")
plt.savefig(f"{new_folder}/predictions")

# Reconstruction
kept_triplets = []
for triplets, prediction in zip(test_triplets_complete, predictions[:, 1]):
    if prediction > 0.1:
        kept_triplets.append(triplets)

reco_xplets_simplified_LUXE(kept_triplets, new_folder, fit="chi squared lin track")
# track_reconstruction_efficiency_simplified_LUXE(f"{new_folder}/reco_xplet_list.npy",
#                                                 "../../qubo_ml_datasets/xi_5/e0gpc_5.0_0000_sl_gen_xplet_list.npy")

track_reconstruction_efficiency_simplified_LUXE(f"{new_folder}/reco_xplet_list_ambiguity_solved.npy",
                                                "../../qubo_ml_datasets/xi_5/e0gpc_5.0_0000_sl_gen_xplet_list.npy")
