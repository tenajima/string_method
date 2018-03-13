from string_method import String
import numpy as np
from example_potential_model import ExamplePotentialModel
import matplotlib.pyplot as plt

def plot_string(dir_name, number_of_point, with_contour=True):
    '''
    ストリングの実行結果を登校線と共に表示
    :param dir_name: solute_num_xx.datまでのpathを指定(必須)
    :param number_of_point: ストリングの点の数を指定(必須)
    :param with_contour: 等高線を表示するか(オプション)
    :return:
    '''
    x_array = np.zeros(number_of_point)
    y_array = np.zeros_like(x_array)
    for i in range(number_of_point):
        f=open(f'{dir_name}solute_num_{i:02d}.dat')
        x, y, _ = map(float,f.readline().split())
        x_array[i] = x
        y_array[i] = y
    plt.plot(x_array, y_array)
    if with_contour:
        example = ExamplePotentialModel([[0, 0, 0]])
        example.plot_contour()
    plt.show()


if __name__ == '__main__':
    # number_of_string_points = 11
    # x = np.linspace(-0.5, 0.5, number_of_string_points)
    # y = np.ones_like(x) * 0.5
    # points = []
    # for i in range(number_of_string_points):
    #     points.append(ExamplePotentialModel([[x[i], y[i], 0]]))
    #     print(points[i].solute_array)
    # example_string = String(points=points, work_dir='example_string/')
    # example_string.stringMethod()

    plot_string(dir_name='./example_string/point_num_11/data_1101/', number_of_point=11)