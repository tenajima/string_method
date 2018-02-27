import numpy as np
from point_of_string import PointOfString
import matplotlib.pyplot as plt

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
            x = self.solute_array[0][0]
        if y is None:
            y = self.solute_array[0][1]
        return (1 - x ** 2 - y **2) ** 2 + y ** 2 / (x ** 2 + y ** 2)

    def get_gradient(self, x=None, y=None):
        '''
        自身のsolute_arrayに関するgradientを返す.
        :return:
        '''
        if x is None:
            x = self.solute_array[0][0] + 1e-10
        if y is None:
            y = self.solute_array[0][1] + 1e-10
        circle = x ** 2 + y ** 2
        return np.array([[-2 * x * (-2 * (-1 + circle) + y**2 / circle**2),
                          2 * y * ( 2 * (-1 + circle) + 1 / circle * (1 - y**2 / circle)),
                          0]])

    def time_development(self, dt):
        '''
        最急降下法により1step時間発展させる
        :param dt:
        :return:
        '''
        self.solute_array += - self.get_gradient() * dt

    def export_solute_array(self,file_name):
        f=open(file_name, 'w')
        for i in range(len(self.solute_array)):
            f.write(f'{self.solute_array[i][0]:06f} {self.solute_array[i][1]:06f} {self.solute_array[i][0]:06f} \n')

    def plot_contour(self, show=True):
        x_ = np.arange(-1.2, 1.2, 0.01)
        y_ = np.arange(-0.2, 1.2, 0.01)
        x, y = np.meshgrid(x_, y_)
        z = self.get_potential_value(x=x, y=y)
        plt.contour(x, y, z, levels=np.arange(-1, 10, 0.1))
        if show:
            plt.show()



if __name__ == '__main__':
    example = ExamplePotentialModel([[0.8, 0.8, 0]])
    # for _ in range(10000):
    #     example.time_development(dt=0.01)
    # print(example.solute_array)
    example.plot_contour()