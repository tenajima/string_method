import numpy as np
from point_of_string import PointOfString

class ExamplePotentialModel(PointOfString):
    '''
    参考論文
    Simplified and improved string method for computing the minimum energy paths in barrier-crossing events
    E, Ren, and Vanden-Eijnden
    J. Chem. Phys. 126, 164103  2007
    参考論文のD. Illustrative exampleの点
    '''
    def __init__(self, solute_array):
        '''
        ここではsolute_arrayの形を1*2のリストで指定する
        :param solute_array:
        '''
        if len(solute_array) is not 2:
            raise TypeError
        self.solute_array = np.array(solute_array)

    def get_potential_value(self, x=None, y=None):
        '''
        ポテンシャルの値を返す.
        x,yに任意の数値を入れればその値を引数にしたポテンシャルの値を返す.
        差分法等で確かめるときに有用.
        デフォルトでは自身が持つsolute_arrayに関するポテンシャルエネルギーを返す.
        :param x:
        :param y:
        :return:
        '''
        if x is None:
            x = self.solute_array[0]
        if y is None:
            y = self.solute_array[1]
        return (1 - x ** 2 - y **2) ** 2 + y ** 2 / (x ** 2 + y ** 2)

    def get_gradient(self, x=None, y=None):
        '''
        自身のsolute_arrayに関するgradientを返す.
        :return:
        '''
        if x is None:
            x = self.solute_array[0]
        if y is None:
            y = self.solute_array[1]
        circle = x ** 2 + y ** 2
        return np.array([-2 * x * (-2 * (-1 + circle) +  y**2 / circle**2),
                          2 * y * ( 2 * (-1 + circle) + 1 / circle * (1 - y**2 / circle))])


    def time_development(self, dt):
        '''
        最急降下法により1step時間発展させる
        :param dt:
        :return:
        '''
        self.solute_array += - self.get_gradient() * dt



if __name__ == '__main__':
    example = ExamplePotentialModel([0.8, 0.8])
    for _ in range(10000):
        example.time_development(dt=0.001)
    print(example.solute_array)
