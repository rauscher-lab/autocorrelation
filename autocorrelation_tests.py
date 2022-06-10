import unittest
import numpy as np

from autocorrelation import (
        one_hot_encode,
        index_autocorrelate,
        index_autocorrelate_old,
        vector_autocorrelate,
        vector_autocorrelate_signal
    )


class TestACMethods(unittest.TestCase):


    def test_oh(self):
        zero = 0
        max_num = 5
        num = 4
        
        self.assertTrue(np.all(one_hot_encode(num, max_num) == np.array([0, 0, 0, 1, 0])))
        self.assertTrue(np.all(one_hot_encode(zero, max_num) == np.array([0, 0, 0, 0, 0])))
        self.assertTrue(np.all(one_hot_encode(zero, 1) == np.array([0])))
        self.assertTrue(np.all(one_hot_encode(1, 1) == np.array([1])))
        self.assertTrue(np.all(one_hot_encode(0, 0) == np.array([0])))
        
        
    def test_id_ac(self):
        index_ts = [35, 35, 36, 0]
        self.assertTrue(np.all(np.abs(index_autocorrelate(index_ts) - np.array([3/4, 1/3, 0, 0])) < 0.000001))
        self.assertTrue(np.all(np.abs(index_autocorrelate([35, 35, 36]) - np.array([1, 1/2, 0])) < 0.000001))
        self.assertTrue(np.all(np.abs(index_autocorrelate([0, 0]) - np.array([0, 0])) < 0.000001))
        self.assertTrue(np.all(np.abs(index_autocorrelate([35, 35, 35]) - np.array([1,1,1])) < 0.000001))
        
        
    def test_id_ac_old(self):
        index_ts = [35, 35, 36, 0]
        self.assertTrue(np.all(np.abs(index_autocorrelate_old(index_ts) - np.array([3/4, 1/3, 0, 0])) < 0.000001))

        
    def test_vec_ac(self):
        ts_array = np.array([
            [1,0,0],[1,0,0],[0.8, 0.6, 0],[0,1,0]
        ])
        self.assertTrue(np.all(np.abs(vector_autocorrelate(ts_array) - np.array([1, 0.8, 0.4, 0.0])) < 0.000001))
        
    def test_vec_ac_signal(self):
        ts_array = np.array([
            [1,0,0],[1,0,0],[0.8, 0.6, 0],[0,1,0]
        ])
        self.assertTrue(np.all(np.abs(vector_autocorrelate_signal(ts_array) - np.array([1, 0.8, 0.4, 0.0])) < 0.000001))

    

if __name__ == '__main__':
    unittest.main()
