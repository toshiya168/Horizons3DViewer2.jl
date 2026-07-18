# 9. Gauss--Newton法・10. Levenberg--Marquardt法

# 9. Gauss--Newton法

## 9.1 推定問題

未知状態ベクトルを

$$
\mathbf{x}
=
\begin{bmatrix}
x&y&v_x&v_y
\end{bmatrix}^{T}
$$

とする。

ここで、

- $\mathbf{x}$：推定したい状態ベクトル
- $x,y$：対象物の位置座標
- $v_x,v_y$：対象物の速度成分
- 上付きの $T$：ベクトルまたは行列の転置

である。

観測ベクトルを

$$
\mathbf{z}
$$

とする。

ここで、$\mathbf{z}$ は実際に観測されたTDOAおよびFDOAをまとめたベクトルである。

観測モデルを

$$
\mathbf{z}
=
h(\mathbf{x})
+
\boldsymbol{\epsilon}
$$

とする。

ここで、

- $h(\mathbf{x})$：状態 $\mathbf{x}$ から理想的なTDOA・FDOAを計算する非線形観測関数
- $\boldsymbol{\epsilon}$：観測ノイズベクトル

である。

観測ノイズがなければ、実際の観測値 $\mathbf{z}$ と予測観測 $h(\mathbf{x})$ は一致する。しかし、実際には観測ノイズが含まれるため、両者は完全には一致しない。

そこで、実際の観測値と予測観測との差ができるだけ小さくなるような状態 $\mathbf{x}$ を求める。

---

## 9.2 残差

現在の状態推定値 $\mathbf{x}$ から予測観測

$$
\hat{\mathbf z}
=
h(\mathbf x)
$$

を計算する。

ここで、$\hat{\mathbf z}$ は現在の状態推定値から計算された予測観測ベクトルであり、上付きの記号 $\hat{\phantom{x}}$ は推定値または予測値であることを表す。

実際の観測値と予測観測との差を、残差ベクトル $\mathbf r$ として

$$
\mathbf r
=
\mathbf z-\hat{\mathbf z}
$$

と定義する。

したがって、

$$
\mathbf r
=
\mathbf z-h(\mathbf x)
$$

である。

ここで、

- $\mathbf r$：観測残差ベクトル
- $\mathbf z$：実際の観測ベクトル
- $\hat{\mathbf z}=h(\mathbf x)$：現在の状態推定値から計算した予測観測ベクトル

である。

残差 $\mathbf r$ が小さいほど、現在の状態推定値 $\mathbf x$ が実際の観測をよく説明していることを意味する。

---

## 9.3 最小二乗問題

状態推定では、残差

$$
\mathbf r
=
\mathbf z-h(\mathbf x)
$$

ができるだけ小さくなるような状態 $\mathbf x$ を求める。

単純に残差の各成分を足し合わせると、正の残差と負の残差が打ち消し合う可能性がある。そのため、残差の二乗を用いて大きさを評価する。

観測ノイズの影響を考慮した目的関数を

$$
J(\mathbf x)
=
\frac12
\mathbf r^{T}
R^{-1}
\mathbf r
$$

と定義する。

ここで、

- $J(\mathbf x)$：最小化する目的関数
- $\mathbf r$：観測残差ベクトル
- $R$：観測ノイズの共分散行列
- $R^{-1}$：観測ノイズ共分散行列の逆行列
- $\mathbf r^T$：残差ベクトル $\mathbf r$ の転置

である。

この目的関数は、観測精度を考慮した残差の二乗和であり、**重み付き最小二乗法**と呼ばれる。

### 残差を二乗する理由

残差の各成分を

$$
\mathbf r
=
\begin{bmatrix}
r_1\\
r_2\\
\vdots\\
r_m
\end{bmatrix}
$$

とする。

ここで、

- $m$：観測量の数
- $r_j$：第 $j$ 観測量の残差

である。

観測ノイズが互いに独立で、すべて同じ分散を持つ場合には、

$$
R
=
\sigma^2 I
$$

と書ける。

ここで、

- $\sigma$：各観測ノイズの標準偏差
- $\sigma^2$：各観測ノイズの分散
- $I$：単位行列

である。

このとき、

$$
R^{-1}
=
\frac{1}{\sigma^2}I
$$

であるため、目的関数は

$$
J(\mathbf x)
=
\frac{1}{2\sigma^2}
\mathbf r^T\mathbf r
$$

となる。

さらに、

$$
\mathbf r^T\mathbf r
=
r_1^2+r_2^2+\cdots+r_m^2
$$

であるから、

$$
J(\mathbf x)
=
\frac{1}{2\sigma^2}
\sum_{j=1}^{m}r_j^2
$$

となる。

ここで、$\sum$ は複数の項を足し合わせる総和記号である。

したがって、この目的関数を最小化することは、観測値と予測観測の差の二乗和を最小化することに対応する。

残差を二乗することで、

- 正の残差と負の残差が打ち消し合わない
- 大きな残差をより強く評価できる
- 微分可能な滑らかな目的関数になる

という利点がある。

### なぜ $R^{-1}$ を掛けるのか

観測量ごとに測定精度が異なる場合、すべての残差を同じ重みで扱うことは適切ではない。

例えば、2つの独立な観測残差を

$$
\mathbf r
=
\begin{bmatrix}
r_1\\
r_2
\end{bmatrix}
$$

とし、それぞれの観測ノイズ分散を $\sigma_1^2$、$\sigma_2^2$ とすると、共分散行列は

$$
R
=
\begin{bmatrix}
\sigma_1^2&0\\
0&\sigma_2^2
\end{bmatrix}
$$

となる。

ここで、

- $\sigma_1$：第1観測量のノイズ標準偏差
- $\sigma_2$：第2観測量のノイズ標準偏差

である。

この逆行列は

$$
R^{-1}
=
\begin{bmatrix}
\dfrac{1}{\sigma_1^2}&0\\
0&\dfrac{1}{\sigma_2^2}
\end{bmatrix}
$$

であるため、目的関数は

$$
J(\mathbf x)
=
\frac12
\left(
\frac{r_1^2}{\sigma_1^2}
+
\frac{r_2^2}{\sigma_2^2}
\right)
$$

となる。

したがって、

- 分散が小さい観測、すなわち精度の高い観測には大きな重みが付く
- 分散が大きい観測、すなわち精度の低い観測には小さな重みが付く

ことになる。

つまり、$R^{-1}$ は観測の信頼度を表す重み行列として働く。

TDOAとFDOAでは、観測量の単位も大きさも異なる。TDOAは時間、FDOAは周波数であるため、残差をそのまま二乗して加えると、一方の観測量だけが目的関数を支配する可能性がある。

$R^{-1}$ を用いることで、それぞれの観測誤差を標準偏差で正規化した形で比較できる。

### 観測誤差に相関がある場合

観測誤差に相関がある場合、共分散行列 $R$ は対角行列とはならない。

例えば、2つのTDOA観測誤差に共分散がある場合、

$$
R_{\rm TDOA}
=
\sigma_t^2
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

となる。

この場合、目的関数

$$
\mathbf r^TR^{-1}\mathbf r
$$

は、各残差の大きさだけでなく、残差間の相関も考慮して評価する。

したがって、相関のある2つの観測を、独立な2つの観測として二重に評価することを防ぐことができる。

### なぜ係数 $\dfrac12$ を付けるのか

目的関数の先頭にある

$$
\frac12
$$

は、目的関数を微分したときに現れる係数2を打ち消すために付けられている。

例えば、

$$
J(r)
=
\frac12r^2
$$

を $r$ で微分すると、

$$
\frac{dJ}{dr}
=
r
$$

となる。

一方、係数 $\frac12$ がない場合には、

$$
\frac{d}{dr}r^2
=
2r
$$

となる。

係数 $\frac12$ は計算式を簡潔にするためのものであり、目的関数を最小にする状態 $\mathbf x$ 自体には影響しない。

### 正規分布との関係

観測ノイズが

$$
\boldsymbol{\epsilon}
\sim
\mathcal N(\mathbf 0,R)
$$

に従うと仮定する。

ここで、

- $\mathcal N(\mathbf 0,R)$：平均ベクトルが $\mathbf 0$、共分散行列が $R$ の多変量正規分布
- $\mathbf 0$：すべての成分が0のベクトル

である。

このとき、状態 $\mathbf x$ が与えられた場合に観測値 $\mathbf z$ が得られる確率密度は、

$$
p(\mathbf z\mid\mathbf x)
\propto
\exp
\left[
-\frac12
\mathbf r^TR^{-1}\mathbf r
\right]
$$

と表される。

ここで、

- $p(\mathbf z\mid\mathbf x)$：状態 $\mathbf x$ のもとで観測 $\mathbf z$ が得られる条件付き確率密度
- $\propto$：両辺が比例することを表す
- $\exp(\cdot)$：指数関数

である。

確率密度 $p(\mathbf z\mid\mathbf x)$ を最大にすることは、

$$
\frac12
\mathbf r^TR^{-1}\mathbf r
$$

を最小にすることと等しい。

したがって、観測ノイズが正規分布に従う場合、重み付き最小二乗法によって得られる状態推定値は、観測が得られる確率を最大にする**最尤推定値**に対応する。

最終的に、推定状態は

$$
\boxed{
\hat{\mathbf x}
=
\mathop{\operatorname{argmin}}_{\mathbf x}
J(\mathbf x)
}
$$

として求められる。

ここで、

- $\hat{\mathbf x}$：推定された状態ベクトル
- $\operatorname{argmin}$：目的関数を最小にする変数の値を表す記号

である。

---

## 9.4 線形化

観測関数 $h(\mathbf x)$ は、距離や視線単位ベクトルを含む非線形関数である。そのため、目的関数を直接最小化することは一般に容易ではない。

そこで、現在の状態推定値 $\mathbf x_k$ の近傍で観測関数を一次近似する。

ここで、添字 $k$ は反復計算の回数を表す。

次の状態を

$$
\mathbf x_k+\Delta\mathbf x
$$

とする。

ここで、$\Delta\mathbf x$ は現在の推定値に加える状態修正量である。

観測関数を一次のTaylor展開によって近似すると、

$$
h(\mathbf x_k+\Delta\mathbf x)
\approx
h(\mathbf x_k)
+
H_k\Delta\mathbf x
$$

となる。

ここで、

$$
H_k
=
\left.
\frac{\partial h}{\partial\mathbf x}
\right|_{\mathbf x=\mathbf x_k}
$$

は、現在の推定値 $\mathbf x_k$ における観測関数のヤコビアンである。

ヤコビアン $H_k$ は、状態の各成分がわずかに変化したとき、各観測量がどの程度変化するかを表す行列である。

現在の推定値における残差を

$$
\mathbf r_k
=
\mathbf z-h(\mathbf x_k)
$$

とする。

修正後の状態における残差は

$$
\mathbf r(\mathbf x_k+\Delta\mathbf x)
=
\mathbf z-h(\mathbf x_k+\Delta\mathbf x)
$$

である。

線形化した観測関数を代入すると、

$$
\mathbf r(\mathbf x_k+\Delta\mathbf x)
\approx
\mathbf z-
\left[
h(\mathbf x_k)
+
H_k\Delta\mathbf x
\right]
$$

となる。

したがって、

$$
\mathbf r(\mathbf x_k+\Delta\mathbf x)
\approx
\mathbf r_k-H_k\Delta\mathbf x
$$

となる。

以降では記号を簡潔にするため、$H_k$ を $H$ と書く。

---

## 9.5 正規方程式

線形化された残差

$$
\mathbf r_k-H\Delta\mathbf x
$$

を用いると、状態修正量 $\Delta\mathbf x$ に関する目的関数は

$$
J(\Delta\mathbf x)
=
\frac12
\left(
\mathbf r_k-H\Delta\mathbf x
\right)^T
R^{-1}
\left(
\mathbf r_k-H\Delta\mathbf x
\right)
$$

となる。

この目的関数を $\Delta\mathbf x$ について最小化する。

勾配を0とすると、

$$
-H^TR^{-1}
\left(
\mathbf r_k-H\Delta\mathbf x
\right)
=
\mathbf 0
$$

となる。

これを展開すると、

$$
-H^TR^{-1}\mathbf r_k
+
H^TR^{-1}H\Delta\mathbf x
=
\mathbf 0
$$

となる。

したがって、

$$
\boxed{
H^{T}R^{-1}H
\Delta\mathbf x
=
H^{T}R^{-1}\mathbf r_k
}
$$

を得る。

この式を正規方程式と呼ぶ。

ここで、

- $H^TR^{-1}H$：状態修正量に対する係数行列
- $H^TR^{-1}\mathbf r_k$：残差から計算される右辺ベクトル
- $\mathbf 0$：すべての成分が0のベクトル

である。

$H^{T}R^{-1}H$ が逆行列を持つ場合、

$$
\boxed{
\Delta\mathbf x
=
(H^{T}R^{-1}H)^{-1}
H^{T}R^{-1}\mathbf r_k
}
$$

を得る。

ただし、数値計算では逆行列を明示的に計算するよりも、

$$
H^{T}R^{-1}H
\Delta\mathbf x
=
H^{T}R^{-1}\mathbf r_k
$$

という連立一次方程式を直接解く方が、一般に数値的に安定である。

状態の更新式は

$$
\boxed{
\mathbf x_{k+1}
=
\mathbf x_k+\Delta\mathbf x
}
$$

となる。

この計算を、残差または状態修正量が十分小さくなるまで繰り返す。

---

# 10. Levenberg--Marquardt法
レーベンバーグ・マーカート法と読む。略して LM法（エルエム法） と呼ばれることも多い。

## 10.1 Gauss--Newton法の問題点

Gauss--Newton法は、現在の推定値の近傍で観測関数を線形化し、その近似に基づいて状態を更新する。

現在の推定値が真値に十分近い場合には、高速に収束することが期待できる。

一方で、

- 初期値が真値から大きく離れている
- $H^{T}R^{-1}H$ が悪条件である
- 観測関数の非線形性が強い
- 受信局の幾何配置が悪く、状態を十分に識別できない

場合には、更新量が過大になったり、正しくない方向に更新されたりして、発散することがある。

ここで、悪条件とは、入力のわずかな変化によって解が大きく変化しやすい状態を意味する。

---

## 10.2 ダンピング

Gauss--Newton法を安定化するため、係数行列の対角成分に正の値を加える。

Levenberg--Marquardt法の正規方程式は

$$
\boxed{
\left(
H^{T}R^{-1}H+\lambda I
\right)
\Delta\mathbf x
=
H^{T}R^{-1}\mathbf r_k
}
$$

となる。

ここで、

- $\lambda$：ダンピング係数
- $I$：状態ベクトルの次元に対応する単位行列
- $\lambda I$：更新量が過大になることを抑えるダンピング項

である。

$H^{T}R^{-1}H+\lambda I$ が逆行列を持つ場合、更新量は

$$
\boxed{
\Delta\mathbf x
=
\left(
H^{T}R^{-1}H+\lambda I
\right)^{-1}
H^{T}R^{-1}\mathbf r_k
}
$$

となる。

ただし、実際の数値計算では、Gauss--Newton法と同様に逆行列を明示的に計算せず、連立一次方程式を直接解くことが望ましい。

---

## 10.3 ダンピング係数の役割

### $\lambda$ が小さい場合

$\lambda$ が十分小さい場合、

$$
H^{T}R^{-1}H+\lambda I
\approx
H^{T}R^{-1}H
$$

となる。

そのため、Levenberg--Marquardt法の更新はGauss--Newton法の更新に近づく。

現在の推定値が真値に近く、線形近似が十分に正確な場合には、$\lambda$ を小さくすることでGauss--Newton法の高速な収束性を利用できる。

### $\lambda$ が大きい場合

$\lambda$ が十分大きい場合、

$$
H^{T}R^{-1}H+\lambda I
\approx
\lambda I
$$

となる。

このとき、

$$
\lambda I\Delta\mathbf x
\approx
H^{T}R^{-1}\mathbf r_k
$$

であるから、

$$
\Delta\mathbf x
\approx
\frac1\lambda
H^{T}R^{-1}\mathbf r_k
$$

となる。

目的関数の勾配は

$$
\nabla J(\mathbf x_k)
=
-H^{T}R^{-1}\mathbf r_k
$$

である。

ここで、$\nabla J$ は目的関数 $J$ の勾配ベクトルを表す。

したがって、

$$
\Delta\mathbf x
\approx
-\frac1\lambda
\nabla J(\mathbf x_k)
$$

となり、目的関数が減少する方向である負の勾配方向へ、小さな更新を行うことになる。

このため、$\lambda$ が大きい場合には、Levenberg--Marquardt法は勾配降下法に近い安全な更新となる。

---

## 10.4 本コードでの更新

本コードでは、反復回数 $k$ ごとに次の処理を行う。

1. 現在の状態推定値 $\mathbf x_k$ から予測観測 $h(\mathbf x_k)$ を計算する
2. 残差 $\mathbf r_k=\mathbf z-h(\mathbf x_k)$ を計算する
3. 数値微分によってヤコビアン $H_k$ を計算する
4. LM法の正規方程式を解いて更新量 $\Delta\mathbf x$ を計算する
5. $\mathbf x_{k+1}=\mathbf x_k+\Delta\mathbf x$ によって状態を更新する

ここで、数値ヤコビアンとは、状態の各成分をわずかに変化させたときの観測関数の変化から、偏微分を近似したヤコビアンである。

収束判定には、例えば

$$
\|\Delta\mathbf x\|
$$

または

$$
\|\mathbf r_k\|
$$

が十分小さいことを利用する。

ここで、$\|\cdot\|$ はベクトルのユークリッドノルムを表す。

状態修正量

$$
\|\Delta\mathbf x\|
$$

が十分小さい場合、反復によって状態推定値がほとんど変化しなくなったことを意味する。

残差

$$
\|\mathbf r_k\|
$$

が十分小さい場合、予測観測が実際の観測値によく一致していることを意味する。

---

# まとめ

Gauss--Newton法は、非線形観測関数を現在の推定値の近傍で線形化し、

$$
H^{T}R^{-1}H
\Delta\mathbf x
=
H^{T}R^{-1}\mathbf r_k
$$

を解くことで状態修正量を求める。

$H^{T}R^{-1}H$ が逆行列を持つ場合には、

$$
\Delta\mathbf x
=
(H^{T}R^{-1}H)^{-1}
H^{T}R^{-1}\mathbf r_k
$$

と表せる。

真値に近い適切な初期値が与えられた場合、Gauss--Newton法は高速に収束する。

一方、Levenberg--Marquardt法は、

$$
\left(
H^{T}R^{-1}H+\lambda I
\right)
\Delta\mathbf x
=
H^{T}R^{-1}\mathbf r_k
$$

のようにダンピング項 $\lambda I$ を加える。

これにより、

- $\lambda$ が小さい場合にはGauss--Newton法に近い高速な更新
- $\lambda$ が大きい場合には勾配降下法に近い安全な更新

を行うことができる。

本サンプルコードでは、このLevenberg--Marquardt法を用いて、TDOAとFDOAの観測値から対象物の位置および速度を同時に推定している。