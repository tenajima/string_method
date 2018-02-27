from abc import ABCMeta, abstractmethod
import numpy as np
class PointOfString(metaclass=ABCMeta):

    @abstractmethod
    def __init__(self, solute_array):
        '''
        :param solute_array:
        ストリングの各点における粒子の構成を表す.
        N個の粒子で構成される点の場合,N*3のリストとなる(3はx,y,z座標の次元).
        '''
        raise NotImplementedError


    @abstractmethod
    def time_development(self, dt):
        '''
        ストリング法において各点をどのように時間発展させるかを実装する.
        :param dt:
        :return:
        '''
        raise NotImplementedError()

    @abstractmethod
    def export_solute_array(self):
        '''
        ログを残すための関数を実装する.
        :return:
        '''
        raise NotImplementedError




