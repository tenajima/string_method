# ストリング法
参考論文
Simplified and improved string method for computing the minimum energy paths in barrier-crossing events
E, Ren, and Vanden-Eijnden
J. Chem. Phys. 126, 164103  2007
参考論文のD. Illustrative exampleの点

参考論文のストリング法を実装したもの.

## 必要なライブラリ
- numpy(必須)
- scipy(必須)
- matplotlib(必須)
- OpenMM(ADPのサンプルを実行する際は必須)

## 使い方
まずはストリングの各点を表現するクラスを実装する.
point_of_string.pyのPointOfStringクラスを継承し,継承先でabstractmethodを実装する.実装しないとNotImplementErrorを出す.
そのストリングの各点を表現するクラスを定義した後,リストに各点に関するインスタンスを入れたものを用いてStringクラスをイニシャライズする.Stringクラスのstring_method()を実行すればストリング法が行われる.

## 例
