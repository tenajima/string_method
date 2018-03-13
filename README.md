# ストリング法
参考論文
Simplified and improved string method for computing the minimum energy paths in barrier-crossing events
E, Ren, and Vanden-Eijnden
J. Chem. Phys. 126, 164103  2007

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
### 論文内のD. Illustrative exampleの点の実装について
- example_string.pyにおいてPointOfStringクラスを継承したクラス「ExamplePotentialModel」を定義する.PointOfStringクラスにおけるabstractmethodを継承先でも実装する必要があることに注意する.
- example_string_method.pyにおいてString.pyで定義されているStringクラスを用いてストリング法を実行できる.イニシャライズにはPointOfStringクラスを継承したExamplePotentialModelのインスタンスが格納されているリストによってイニシャライズする.イニシャライズする際にログファイルを残すディレクトリを「work_dir=」で指定できる.
- StringクラスのインスタンスのstringMethod()を実行することでストリング法を行うことができる.

### ADPについて
基本的な手順は同じであるが,ADPを用いる際に注意する点を記述する.
- ADPはストリング法を用いる際の時間刻みdtに敏感なので適切な値を設定してあげる必要がある.
- OpenMMを用いる際には物理量がQuantityクラスによって表されているので単位が文字列でついている.そのためnumpyで計算しようとしてもエラーが出る.この部分を吸収するためにADPクラスが定義されているADPクラスでは,物理量を表す部分にsetterとgetterを定義することで,Stringクラスからは特別な処理をしなくても実装できるようにした.
