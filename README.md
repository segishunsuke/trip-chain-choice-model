このリポジトリでは，トリップチェイン全体の効用を最大化する旅行者の意思決定を明示的に表現した，観光周遊行動モデルを扱うためのプログラムを公開しています．トリップチェインの観測データを用いたモデルの推定や，推定結果を利用したトリップチェインの予測を行うことができます．

プログラムについての問い合わせは，下記のメールアドレスにお願いします．

<img src="./assets/images/mail.png" width="200px">


## Cythonライブラリのインストール

このプログラムを用いるには，Python環境およびCythonライブラリが必要です．

お使いのPython環境にCythonライブラリがインストールされていない場合は，プロンプト上で以下のコマンドを入力してインストールを行って下さい．

```
pip install Cython
```

### Windows環境の注意点

Windows環境でCythonプログラムをコンパイルするには，C/C++のコンパイラ（Visual Studio）の導入が必要です．コンパイラは無償で利用可能です．

以下の手順に従って，必要な機能をインストールして下さい．

1. Visual Studio Communityのインストール

- [Microsoftの公式サイト](https://visualstudio.microsoft.com/downloads/)にアクセスし，「Visual Studio Community」をダウンロードします．
- インストーラを起動し，「Python開発」を選択してインストールします．

2. Build Tools for Visual Studioのインストール

- [同じページ](https://visualstudio.microsoft.com/downloads/)の「Tools for Visual Studio」セクションから「Build Tools for Visual Studio」をダウンロードします．
- インストーラを起動し，「C++によるデスクトップ開発」を選択してインストールします．

## プログラムのダウンロードとコンパイル

[codesフォルダ](./codes)内のファイルを全て同一のフォルダにダウンロードして下さい．ファイルの内容は以下の通りです．

- `trip_chain_simulator.pyx`, `trip_chain_simulator.pxd`: 観光周遊行動モデルを扱うライブラリの本体です．
- `geneticr.pyx`, `geneticr.pxd`: 実数値の遺伝的アルゴリズムによる関数最適化を行うライブラリです．
- `mt19937ar.c`: メルセンヌツイスタによる乱数生成を行うCコードです．[開発者により公開されているコード](https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/mt19937ar.html)を編集したものです．
- `setup.py`: コンパイルに利用するPythonファイルです．
- `execute.py`: `trip_chain_simulator`ライブラリの利用例が載っているPythonファイルです．

ファイルを全てダウンロードしたら，以下のコマンドを実行してプログラムをコンパイルして下さい．

```
python setup.py build_ext --inplace
```

コンパイルは一度実行すれば，それ以降は行う必要はありません．


## 設定ファイルとデータファイルの準備

`trip_chain_simulator`ライブラリを使うには，以下の四つの設定ファイルとデータファイルを準備する必要があります．これらのファイルは全てカンマ区切りのCSVファイルとして準備して下さい．
テストデータ用のファイルは[test-dataフォルダ](./test-data)に置かれています．

### 設定ファイル

設定ファイルは以下のような4行2列構成のCSVファイルとして準備して下さい．
[test-dataフォルダ](./test-data)内の`input settings.csv`が例です．
1列目は項目名，2列目は設定値です．
項目の順序は変更できません．指定された順序を守って下さい．

| Number of places | 11 |
| :---: | :---: |
| Number of ports | 2 |
| Shift parameter of Poisson likelihood | 0.0001 |
| OD cost normalization | 0 |

各項目の意味は以下の通りです．

- `Number of places`: 旅行者が訪問可能な場所（観光地）の数です．これらの場所はトリップチェインの起点・終点になることはできません．
- `Number of ports`: トリップチェインの起点・終点となることができる場所（空港など，以下では便宜上「空港」と呼称）の数です．これらの場所は旅行者の訪問の対象となることはできません．
- `Shift parameter of Poisson likelihood`: ポアソン疑似尤度を評価する際に，ゼロ予測値を避けるための補正項．通常は0.0001に設定することを推奨します．
- `OD cost normalization`: 旅行費用の単位が分析に影響を与えることを防ぐため，OD間旅行費用はこの設定値で除されて基準化されたうえで，プログラムに利用されます．この設定値にゼロを指定すると，OD間旅行費用の95%タイル値が基準化に使われます．

Number of placesを$`K`$，Number of portsを$`M`$とするとき，各観光地および空港には以下のIDが割り当てられます．

- 観光地: 0, 1, 2, ..., $`K - 1`$
- 空港: $`K`$, $`K + 1`$, ..., $`K + M - 1`$

### OD間旅行費用のデータファイル

OD間旅行費用のデータファイルは，観光地・空港間の旅行費用を記録した，以下のような3列構成のヘッダ付きCSVファイルとして準備して下さい．
[test-dataフォルダ](./test-data)内の`input od cost.csv`が例です．
1行目はヘッダです．
1列目は起点のID，2列目は終点のID，3列目がOD間の旅行費用です．
列の順序は変更できません．指定された順序を守って下さい．

| origin | destination | cost |
| :---: | :---: | :---: |
| 0 | 1 | 1.0 |
| 0 | 2 | 2.0 |
| 0 | 3 | 3.0 |
| $`\vdots`$ | $`\vdots`$ | $`\vdots`$ |

特定のODペアがファイルに含まれていない場合，その旅行費用はゼロとして扱われます．
異なる空港間の旅行費用は設定可能ですが，プログラム中で使われることはありません．

### 観測されたトリップチェインのデータファイル

観測されたトリップチェインのデータファイルは，以下のようなヘッダ付きCSVファイルとして準備して下さい．
列数は任意であり，トリップチェインのデータは何列目にあっても構いません．

| $`\cdots`$ | trip chain | $`\cdots`$ |
| :---: | :---: | :---: |
| $`\cdots`$ | \[10 4 10\] | $`\cdots`$ |
| $`\cdots`$ | \[10 3 6 10\] | $`\cdots`$ |
| $`\cdots`$ | \[10 2 4 7 11\] | $`\cdots`$ |
| $`\vdots`$ | $`\vdots`$ | $`\vdots`$ |

トリップチェインは，訪問した観光地・空港のIDを半角スペースで区切り，全体を\[ \]で囲んで表現します．
トリップチェインの先頭および末尾のID（起点と終点）は，空港のIDでなければいけません．先頭と末尾の他のIDは，観光地のIDでなければいけません．
また，同一の観光地を複数回訪問する経路は，このプログラムではサポートされていません．

### パラメータ値の設定ファイル

パラメータ値の設定ファイルは，以下のような2列構成のCSVファイルとして準備して下さい．
行数はNumber of placesを$`K`$としたとき，$`K + 2`$行になります．
[test-dataフォルダ](./test-data)内の`input initial parameter values.csv`が例です．
1列目は項目名，2列目は設定値です．
項目の順序は変更できません．指定された順序を守って下さい．

| alpha\[0\] | -1.0 |
| :---: | :---: |
| alpha\[1\] | 0.0 |
| $`\vdots`$ | $`\vdots`$ |
| alpha\[K-1\] | -1.0 |
| beta | 7.0 |
| sigma_t | 0.3 |

各パラメータの意味を説明するために，このプログラムにおける旅行者のトリップチェイン選択モデルを紹介します．
このプログラムでは，個々の旅行者は以下の効用最大化問題の解をヒューリスティックにより求めてトリップチェインを選択すると仮定しています．

```math
\max_{n, x_1, x_2, ..., x_n } \sum_{j=1}^{n} \left( \alpha_{x_j} + \varepsilon_{x_j} \right) - \beta \sum_{j=0}^n \bar{c}_{x_j, x_{j+1}} \left( 1 + \xi_{x_j, x_{j+1}} \right)
```
```math
\mathrm{s.t.} \: n \ge 1
```
```math
x_0 = \mathrm{port \: of \: entrance}, x_{n+1} = \mathrm{port \: of \: departure}
```
```math
x_j \in \{ 0, \cdots, K - 1 \} \quad (1 \le j \le n)
```
```math
x_j \neq x_{r} \quad (1 \le j \le n, \: 1 \le r \le n, \: j \neq r)
```

ここで，$`n`$は訪問する観光地の数，$`x_j`$は$`j`$番目に訪問する観光地のID，$`\alpha_{k} + \varepsilon_{k}`$は観光地$`k`$の訪問効用，$`\bar{c}_{k, l} ( 1 + \xi_{k, l} )`$は$`kl`$間の旅行費用です．
$`\alpha_{k}`$は観光地$`k`$の訪問効用の期待値であり，パラメータです．
$`\beta`$は旅行費用を効用の単位に変換する係数であり，パラメータです．
$`\bar{c}_{k, l}`$は，OD間旅行費用のデータファイルに記載された$`kl`$間の旅行費用を，設定ファイルの`OD cost normalization`で除した数値です．
$`\varepsilon_{k}`$と$`\xi_{k, l}`$は，平均が0の独立な正規分布に従う確率変数です．正規分布の標準偏差は，それぞれ1と$`\sigma_t`$です．$`\sigma_t`$はパラメータです．
効用最大化問題を解く時点では，$`\varepsilon_{k}`$と$`\xi_{k, l}`$の実現値は確定しています．

以上を踏まえ，パラメータ値の設定ファイルの項目は以下の通りです．

- `alpha[k]`: IDがkの観光地の訪問効用の期待値$`\alpha_{k}`$
- `beta`: 旅行費用抵抗係数$`\beta`$
- `sigma_t`: 旅行費用のばらつきの大きさ$`\sigma_t`$

パラメータの推定を行う場合も設定ファイルは必須です．最適化の初期値として，記載されたパラメータ値が利用されます．

## プログラムの実行

プログラムを使用する際は、`trip_chain_simulator.pyx`の置かれているフォルダをカレントディレクトリにして下さい。
プログラムの実行方法は以下の通りです。

### オブジェクトの作成とインプットファイルの読み込み

まず、`trip_chain_simulator`をインポートします。

```
import trip_chain_simulator
```

次に、`trip_chain_simulator`ライブラリの`Trip_chain_simulator`クラスのオブジェクトを作成します。
```
simulator = trip_chain_simulator.Trip_chain_simulator()
```

次に、準備した設定ファイル・データファイルをオブジェクトに読み込ませます。
```
simulator.read_input_data("input settings.csv", "input od cost.csv", "input trip chain data.csv", "trip chain", "input initial parameter values.csv")
```
このメソッドの引数は、順番に、設定ファイル名、OD間旅行費用のデータファイル名、

### パラメータの推定



