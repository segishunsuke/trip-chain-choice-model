このリポジトリでは，トリップチェイン全体の効用を最大化する旅行者の意思決定を明示的に表現した，観光周遊行動モデルを扱うためのプログラムを公開しています。トリップチェインの観測データを用いたモデルの推定や、推定結果を利用したトリップチェインの予測を行うことができます。

プログラムについての問い合わせは，下記のメールアドレスにお願いします．

<img src="./assets/images/mail.png" width="200px">


## Cythonライブラリのインストール

このプログラムを用いるには，Python環境およびCythonライブラリが必要です。

お使いのPython環境にCythonライブラリがインストールされていない場合は，プロンプト上で以下のコマンドを入力してインストールを行って下さい。

```
pip install Cython
```

### Windows環境の注意点

Windows環境でCythonプログラムをコンパイルするには，C/C++のコンパイラ（Visual Studio）の導入が必要です。コンパイラは無償で利用可能です。

以下の手順に従って，必要な機能をインストールしてください。

1. Visual Studio Communityのインストール

- [Microsoftの公式サイト](https://visualstudio.microsoft.com/downloads/)にアクセスし、「Visual Studio Community」をダウンロードします。
- インストーラを起動し，「Python開発」を選択してインストールします。

2. Build Tools for Visual Studioのインストール

- [同じページ](https://visualstudio.microsoft.com/downloads/)の「Tools for Visual Studio」セクションから「Build Tools for Visual Studio」をダウンロードします。
- インストーラを起動し，「C++によるデスクトップ開発」を選択してインストールします。





このプログラムが用いる，河床標高の設定方法を述べます．このプログラムは，広矩形単断面を持つ開水路の不等流計算の基礎式である，
```math
\frac{dH}{dx} + \frac{1}{2g} \frac{d}{dx} \left( \frac{Q}{Bh} \right)^2 + \frac{n^2 Q^2}{B^2 h^{10/3}} = 0
```
を用いています．ここで，$`H`$は水面標高(m)を