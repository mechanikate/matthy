import pytest
from mattrix import *
class TestMatrix:
    def test_dim_0x0(self): # make sure 0x0 matrices dimensions return [0,0]
        assert dim(M([[]])) == [0,0]
        assert dim(M([[],[],[]])) == [0,0]
    def test_det_size(self): # make sure non-square matrices don't work with M.det
        with pytest.raises(ValueError):
            M([[1,2,3],[6,5,4]]).det()
    def test_det_2x2(self): # check 2x2 matrices work with M.det
        assert M([[16,24],[2,1]]).det() == -32
    def test_det_bigger(self): # check M.det works with bigger matrices
        assert M([[5,2,19,22],[49,102,194,4],[432,9,3,930],[23,2039,39,203]]).det() == -3691612155
    def test_mul_1x1(self): # 1x1 mul test
        assert M([[15]]).mul(M([[6]])).v == [[90]]
    def test_mul_bigger(self): # random mul test
        assert M([[19,321,453],[543,19,192]]).mul(M([[432,49],[109,194],[555,875]])).v == [[294612,459580],[343207,198293]]
    def test_mul_size(self): # make sure multiplication doesn't multiply incompatible matrices:!
        with pytest.raises(ValueError):
            M([[5,2]]).mul(M([[2,5]]))
    def test_mul_eye(self):
        for m in range(1,5): # test multiplicative identity working correctly
            for n in range(1,5):
                assert eye(m,m).mul(ones(m,n)).v == ones(m,n).v
    def test_add_0x0(self): # adding 2 0x0 matrices should yield a 0x0 matrix
        assert M([[]]).add(M([[]])).v == [[]]
    def test_add_bigger(self): # random add test
        assert M([[1,2,3],[4324,54365,643]]).add(M([[654,6535,430940],[4854,65486,546]])).v == [[655,6537,430943],[9178,119851,1189]]
    def test_kmul_0x0(self): # 0x0 scalar multiplication test
        assert M([[]]).kmul(15).v == [[]]
    def test_kmul_bigger(self): # bigger scalar multiplication test
        assert M([[0,569,0],[194,349,549]]).kmul(3).v == [[0,1707,0],[582,1047,1647]]
    def test_sum_0x0(self): # test summing 0x0 matrix
        assert M([[]]).sum() == 0
    def test_sum_bigger(self): # test summing bigger matrices
        assert M([[1,2,3],[4,5,6]]).sum() == 21
        assert M([[1,4],[2,5],[3,6]]).sum() == 21
    def test_min_0x0(self): # test min of 0x0 matrix raises ValueError
        with pytest.raises(ValueError):
            M([[]]).min()
    def test_max_0x0(self): # test max of 0x0 matrix raises ValueError
        with pytest.raises(ValueError):
            M([[]]).max()
    def test_inv_size(self):
        with pytest.raises(ValueError):
            M([[1,2]]).inv()
    def test_inv_1x1(self):
        assert M([[500]]).inv().v == [[0.002]]
    def test_inv_2x2(self):
        assert M([[52,26],[125,2]]).inv().v == [[-1/1573, 1/121], [125/3146, -2/121]]
    def test_pr_linear(self):
        for degree in range(3,6): # should always return in form of y=x (coeff. of 1)
            assert pr([[1,1],[2,2],[3,3],[4,4],[5,5],[0,0]],degree) == [0,1]+([0]*(degree-2))

