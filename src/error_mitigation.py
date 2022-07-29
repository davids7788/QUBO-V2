import qiskit

from qiskit.utils import mitigation
import numpy as np


class ErrorMitigation:
    def __init__(self,
                 backend):
        """Class for creating and Error Mitigation object.
        :param backend: backed object, simulator or Fake Device"""
        self.backend = backend
        self.meas_filter_matrix = None

    def algebraic_mitigation(self, qc):
        """Calibrates the meas filter with the backend and the quantum circuit"""
        cal_circuits, state_labels = mitigation.complete_meas_cal(qr=qc.qregs[0],
                                                                  circlabel='measurement_calibration')
        cal_job = qiskit.execute(cal_circuits,
                                 backend=self.backend,
                                 shots=1024)
        cal_results = cal_job.result()
        meas_fitter = mitigation.CompleteMeasFitter(cal_results,
                                                    state_labels,
                                                    circlabel='measurement_calibration')
        self.meas_filter_matrix = np.linalg.inv(np.array(meas_fitter.cal_matrix))

    @staticmethod
    def dict_to_vector(dictionary):
        """Takes a vqe result eigenstate in dictionary form and transforms it to a vector.
        :param dictionary: dictionary with strings in binary format as keys and probabilities as values
        :return
            vector representing the same information as the dictionary
        """
        len_vector = len(list(dictionary.keys())[0])
        vector = np.zeros(2**len_vector)
        for key, value in dictionary.items():
            vector[int(key, 2)] = value
        return vector

    @staticmethod
    def vector_to_dict(vector):
        """Takes a vector representing states and probabilities and returns a dictionary in the form of vqe eigenstate.
        :param vector: representing eigenstates of vqe result
        :return
            dictionary with strings in binary format as keys and probabilities as values
        """
        def return_binary_string(decimal, length):
            return format(decimal, "b").zfill(length)

        str_length = 0
        for i in range(10):
            if len(vector) == 2**i:
                str_length = i

        dictionary = {}
        for i, value in enumerate(vector):
            if value != 0:
                dictionary.update({return_binary_string(i, str_length): value})
        return dictionary
