import numpy as np
import time
import matplotlib.pyplot as plt
from string_method import String

if __name__ == '__main__':
    print('start string method')
    from ADP import ADP
    start_time = time.time()
    points = []
    number_of_points = 128
    dir_name = f'./vacuum_{number_of_points}/'
    print('pdbファイルを読み込み中')
    count = 0
    for i in range(number_of_points):
        print(f'initial start: {i}')
        point = ADP(number_of_point=i, dir_name=dir_name)
        points.append(point)
    print('count is :', count)
    string = String(points, work_dir=dir_name)
    string.stringMethod()
    elapsed_time = time.time() - start_time
    print('elapsed time is {} [s]'.format(elapsed_time))