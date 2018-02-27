from string_method import String
import numpy as np
from example_potential_model import ExamplePotentialModel
import matplotlib.pyplot as plt


class ExampleString(String):
    def interpolate(self):
        s = [0.0]
        alpha = []
        for i in range(1, self.number_of_points):
            length = 0.0
            for solute_num in range(self.solute_sum):
                for coordinate in range(3):
                    length += (self.Points[i].solute_array[solute_num][coordinate] -
                               self.Points[i - 1].solute_array[solute_num][coordinate]) ** 2.0
            s.append(length ** 0.50 + s[i - 1])
        for i in range(self.number_of_points):
            alpha.append(s[i] / s[self.number_of_points - 1])

        for solute_num in range(self.solute_sum):
            for coordinate in range(3):
                measure = []
                for i in range(self.number_of_points):
                    measure.append(self.Points[i].solute_array[solute_num][coordinate])
                cubicInterp = interp1d(alpha, measure, kind='cubic')

                for i in range(self.number_of_points):
                    self.Points[i].solute_array[solute_num][coordinate] = cubicInterp(i / (self.number_of_points - 1))

def plot_string(dir_name, number_of_point):
    x_array = np.zeros(number_of_point)
    y_array = np.zeros_like(x_array)
    for i in range(number_of_point):
        f=open(f'{dir_name}solute_num_{i:02d}.dat')
        x, y, _ = map(float,f.readline().split())
        x_array[i] = x
        y_array[i] = y
    plt.plot(x_array, y_array)
    example = ExamplePotentialModel([[0, 0, 0]])
    example.plot_contour(show=False)
    plt.show()


if __name__ == '__main__':
    number_of_string_points = 11
    x = np.linspace(-0.5, 0.5, number_of_string_points)
    y = np.ones_like(x) * 0.5
    points = []
    for i in range(number_of_string_points):
        points.append(ExamplePotentialModel([[x[i], y[i], 0]]))
        print(points[i].solute_array)
    example_string = String(points=points, work_dir='example_string/')
    example_string.stringMethod()

    # plot_string(dir_name='./example_string/point_num_11/data_1101/', number_of_point=11)