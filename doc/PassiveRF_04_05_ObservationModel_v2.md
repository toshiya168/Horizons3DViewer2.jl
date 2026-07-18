# 4. 観測モデルと5. 観測ノイズ

## 4. 観測モデル

対象物の状態ベクトルを

$$
\mathbf{x}
=
\begin{bmatrix}
x\\
y\\
v_x\\
v_y
\end{bmatrix}
$$

とする。

ここで、

- $\mathbf{x}$：対象物の状態ベクトル
- $x,y$：対象物の位置座標
- $v_x,v_y$：対象物の速度成分

である。

対象物の位置ベクトル $\mathbf{p}$ と速度ベクトル $\mathbf{v}$ は、それぞれ

$$
\mathbf{p}
=
\begin{bmatrix}
x\\
y
\end{bmatrix},
\qquad
\mathbf{v}
=
\begin{bmatrix}
v_x\\
v_y
\end{bmatrix}
$$

である。

受信局1を基準局とし、TDOAおよびFDOAをまとめた観測ベクトルを

$$
\mathbf{z}
=
\begin{bmatrix}
\Delta t_{21}\\
\Delta t_{31}\\
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

とする。

ここで、

- $\mathbf{z}$：観測ベクトル
- $\Delta t_{21}=t_2-t_1$：受信局2と受信局1の到来時刻差
- $\Delta t_{31}=t_3-t_1$：受信局3と受信局1の到来時刻差
- $\Delta f_{21}=f_2-f_1$：受信局2と受信局1の到来周波数差
- $\Delta f_{31}=f_3-f_1$：受信局3と受信局1の到来周波数差
- $t_i$：受信局 $i$ における受信時刻
- $f_i$：受信局 $i$ における受信周波数

である。

対象物の状態 $\mathbf{x}$ から、理想的な観測量を計算する関数を観測関数 $h$ と呼ぶ。

観測ノイズを考えない理想的な観測モデルは

$$
\mathbf{z}
=
h(\mathbf{x})
$$

と表される。

ここで、

- $h$：状態ベクトルから観測ベクトルを計算する観測関数
- $h(\mathbf{x})$：状態 $\mathbf{x}$ に対応する理想的な観測ベクトル

である。

観測関数 $h$ は、TDOA観測関数とFDOA観測関数から構成される。

---

### 4.1 TDOA観測関数

対象物と受信局 $i$ の間の距離を

$$
r_i
=
\|\mathbf{p}-\mathbf{s}_i\|
$$

とする。

ここで、

- $r_i$：対象物と受信局 $i$ の間の距離
- $\mathbf{p}$：対象物の位置ベクトル
- $\mathbf{s}_i$：受信局 $i$ の位置ベクトル
- $\|\cdot\|$：ベクトルのユークリッドノルム

である。

TDOAに対応する観測関数を $h_{\rm TDOA}$ とすると、

$$
h_{\rm TDOA}(\mathbf{x})
=
\begin{bmatrix}
\dfrac{r_2-r_1}{c}\\
\dfrac{r_3-r_1}{c}
\end{bmatrix}
$$

となる。

ここで、

- $h_{\rm TDOA}$：TDOA観測量を計算する観測関数
- $c$：光速
- $r_2-r_1$：受信局2と受信局1に対する距離差
- $r_3-r_1$：受信局3と受信局1に対する距離差

である。

したがって、

$$
h_{\rm TDOA}(\mathbf{x})
=
\begin{bmatrix}
\Delta t_{21}\\
\Delta t_{31}
\end{bmatrix}
$$

に対応する。

---

### 4.2 FDOA観測関数

受信局 $i$ から対象物へ向かう視線単位ベクトルを

$$
\mathbf{u}_i
=
\frac{\mathbf{p}-\mathbf{s}_i}
{\|\mathbf{p}-\mathbf{s}_i\|}
$$

とする。

ここで、

- $\mathbf{u}_i$：受信局 $i$ から対象物へ向かう視線単位ベクトル
- $\mathbf{p}-\mathbf{s}_i$：受信局 $i$ から対象物へ向かう相対位置ベクトル

である。

受信局が静止している場合、受信局 $i$ に対する距離速度は

$$
\dot r_i
=
\mathbf{u}_i^T\mathbf{v}
$$

である。

ここで、

- $\dot r_i$：対象物と受信局 $i$ の間の距離速度
- $\mathbf{v}$：対象物の速度ベクトル
- 上付きの $T$：ベクトルまたは行列の転置

を表す。

FDOAに対応する観測関数を $h_{\rm FDOA}$ とすると、

$$
h_{\rm FDOA}(\mathbf{x})
=
-\frac{f_0}{c}
\begin{bmatrix}
(\mathbf{u}_2-\mathbf{u}_1)^T\mathbf{v}\\
(\mathbf{u}_3-\mathbf{u}_1)^T\mathbf{v}
\end{bmatrix}
$$

となる。

ここで、

- $h_{\rm FDOA}$：FDOA観測量を計算する観測関数
- $f_0$：対象物が送信する電波の送信周波数
- $\mathbf{u}_2-\mathbf{u}_1$：受信局2と受信局1の視線単位ベクトルの差
- $\mathbf{u}_3-\mathbf{u}_1$：受信局3と受信局1の視線単位ベクトルの差

である。

したがって、

$$
h_{\rm FDOA}(\mathbf{x})
=
\begin{bmatrix}
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

に対応する。

---

### 4.3 全観測関数

TDOAとFDOAをまとめると、全観測関数は

$$
\boxed{
h(\mathbf{x})
=
\begin{bmatrix}
\dfrac{r_2-r_1}{c}\\
\dfrac{r_3-r_1}{c}\\
-\dfrac{f_0}{c}
(\mathbf{u}_2-\mathbf{u}_1)^T\mathbf{v}\\
-\dfrac{f_0}{c}
(\mathbf{u}_3-\mathbf{u}_1)^T\mathbf{v}
\end{bmatrix}
}
$$

となる。

この観測関数は、状態ベクトル

$$
\mathbf{x}
=
\begin{bmatrix}
x\\
y\\
v_x\\
v_y
\end{bmatrix}
$$

から、理想的な観測量

$$
h(\mathbf{x})
=
\begin{bmatrix}
\Delta t_{21}\\
\Delta t_{31}\\
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

を計算する非線形関数である。

観測関数が非線形になるのは、距離

$$
r_i
=
\|\mathbf{p}-\mathbf{s}_i\|
$$

および視線単位ベクトル

$$
\mathbf{u}_i
=
\frac{\mathbf{p}-\mathbf{s}_i}
{\|\mathbf{p}-\mathbf{s}_i\|}
$$

が、対象物の位置 $\mathbf{p}$ に対して非線形に依存するためである。

---

## 5. 観測ノイズ

実際の受信時刻および受信周波数には、受信機雑音、時刻同期誤差、周波数測定誤差などが含まれる。

ここでは、各受信局における時刻測定誤差が互いに独立で、平均0、分散 $\sigma_t^2$ の正規分布に従うと仮定する。

受信局 $i$ の時刻測定誤差を $e_{t,i}$ とすると、

$$
e_{t,i}
\sim
\mathcal{N}(0,\sigma_t^2)
$$

である。

ここで、

- $e_{t,i}$：受信局 $i$ の時刻測定誤差
- $\mathcal{N}(\mu,\sigma^2)$：平均 $\mu$、分散 $\sigma^2$ の正規分布
- $\sigma_t$：各受信局における時刻誤差の標準偏差
- $\sigma_t^2$：各受信局における時刻誤差の分散
- 記号 $\sim$：「その確率分布に従う」こと

を表す。

時刻誤差を含む受信時刻を

$$
\tilde t_i
=
t_i+e_{t,i}
$$

とする。

ここで、

- $\tilde t_i$：時刻誤差を含む受信時刻
- $t_i$：誤差を含まない理想的な受信時刻

である。

受信局2と受信局1のTDOA観測値は

$$
\widetilde{\Delta t}_{21}
=
\tilde t_2-\tilde t_1
$$

である。

ここで、上付きの波線 $\widetilde{\phantom{x}}$ は、観測ノイズを含む測定値であることを表す。

これを展開すると、

$$
\widetilde{\Delta t}_{21}
=
(t_2+e_{t,2})
-
(t_1+e_{t,1})
$$

となり、

$$
\widetilde{\Delta t}_{21}
=
\Delta t_{21}
+
e_{t,2}-e_{t,1}
$$

が得られる。

したがって、TDOAに含まれる観測誤差を

$$
\epsilon_{t,21}
=
e_{t,2}-e_{t,1}
$$

と定義できる。

ここで、$\epsilon_{t,21}$ は受信局2と受信局1のTDOAに含まれる誤差である。

各局の時刻誤差が独立であるため、その分散は

$$
\operatorname{Var}
\left(
\epsilon_{t,21}
\right)
=
\operatorname{Var}
\left(
e_{t,2}-e_{t,1}
\right)
$$

となる。

ここで、$\operatorname{Var}(\cdot)$ は確率変数の分散を表す。

独立な確率変数の分散は加算されるため、

$$
\operatorname{Var}
\left(
\epsilon_{t,21}
\right)
=
\operatorname{Var}(e_{t,2})
+
\operatorname{Var}(e_{t,1})
$$

となり、

$$
\operatorname{Var}
\left(
\widetilde{\Delta t}_{21}
\right)
=
2\sigma_t^2
$$

が得られる。

同様に、

$$
\operatorname{Var}
\left(
\widetilde{\Delta t}_{31}
\right)
=
2\sigma_t^2
$$

となる。

---

### 5.1 TDOA誤差間の相関

2つのTDOA誤差は

$$
\epsilon_{t,21}
=
e_{t,2}-e_{t,1}
$$

および

$$
\epsilon_{t,31}
=
e_{t,3}-e_{t,1}
$$

である。

どちらにも基準局1の時刻誤差 $e_{t,1}$ が含まれているため、2つのTDOA誤差は互いに独立ではない。

両者の共分散は

$$
\operatorname{Cov}
\left(
\epsilon_{t,21},
\epsilon_{t,31}
\right)
=
\sigma_t^2
$$

となる。

ここで、$\operatorname{Cov}(X,Y)$ は確率変数 $X$ と $Y$ の共分散を表す。

したがって、TDOA観測誤差ベクトルを

$$
\boldsymbol{\epsilon}_{\rm TDOA}
=
\begin{bmatrix}
\epsilon_{t,21}\\
\epsilon_{t,31}
\end{bmatrix}
$$

とすると、その共分散行列は

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

ここで、

- $\boldsymbol{\epsilon}_{\rm TDOA}$：TDOA観測誤差ベクトル
- $R_{\rm TDOA}$：TDOA観測誤差の共分散行列

である。

---

### 5.2 FDOA観測誤差

各受信局における周波数測定誤差を $e_{f,i}$ とする。

各局の周波数誤差が互いに独立で、平均0、分散 $\sigma_f^2$ の正規分布に従うと仮定すると、

$$
e_{f,i}
\sim
\mathcal{N}(0,\sigma_f^2)
$$

である。

ここで、

- $e_{f,i}$：受信局 $i$ の周波数測定誤差
- $\sigma_f$：各受信局における周波数誤差の標準偏差
- $\sigma_f^2$：各受信局における周波数誤差の分散

である。

FDOA観測誤差を

$$
\epsilon_{f,21}
=
e_{f,2}-e_{f,1}
$$

および

$$
\epsilon_{f,31}
=
e_{f,3}-e_{f,1}
$$

とすると、TDOAの場合と同様に、

$$
\operatorname{Var}
\left(
\epsilon_{f,21}
\right)
=
2\sigma_f^2
$$

$$
\operatorname{Var}
\left(
\epsilon_{f,31}
\right)
=
2\sigma_f^2
$$

および

$$
\operatorname{Cov}
\left(
\epsilon_{f,21},
\epsilon_{f,31}
\right)
=
\sigma_f^2
$$

となる。

FDOA観測誤差ベクトルを

$$
\boldsymbol{\epsilon}_{\rm FDOA}
=
\begin{bmatrix}
\epsilon_{f,21}\\
\epsilon_{f,31}
\end{bmatrix}
$$

とすると、その共分散行列は

$$
R_{\rm FDOA}
=
\sigma_f^2
\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}
$$

となる。

ここで、

- $\boldsymbol{\epsilon}_{\rm FDOA}$：FDOA観測誤差ベクトル
- $R_{\rm FDOA}$：FDOA観測誤差の共分散行列

である。

---

### 5.3 全観測共分散行列

TDOA観測誤差とFDOA観測誤差が互いに独立であると仮定する。

全観測誤差ベクトルを

$$
\boldsymbol{\epsilon}
=
\begin{bmatrix}
\boldsymbol{\epsilon}_{\rm TDOA}\\
\boldsymbol{\epsilon}_{\rm FDOA}
\end{bmatrix}
=
\begin{bmatrix}
\epsilon_{t,21}\\
\epsilon_{t,31}\\
\epsilon_{f,21}\\
\epsilon_{f,31}
\end{bmatrix}
$$

とする。

ここで、$\boldsymbol{\epsilon}$ はTDOAとFDOAのすべての観測誤差をまとめたベクトルである。

全観測誤差の共分散行列を $R$ とすると、

$$
R
=
\begin{bmatrix}
R_{\rm TDOA} & 0_{2\times2}\\
0_{2\times2} & R_{\rm FDOA}
\end{bmatrix}
$$

となる。

ここで、

- $R$：全観測誤差の共分散行列
- $0_{2\times2}$：すべての要素が0である $2\times2$ 零行列

である。

具体的には、

$$
R
=
\begin{bmatrix}
2\sigma_t^2 & \sigma_t^2 & 0 & 0\\
\sigma_t^2 & 2\sigma_t^2 & 0 & 0\\
0 & 0 & 2\sigma_f^2 & \sigma_f^2\\
0 & 0 & \sigma_f^2 & 2\sigma_f^2
\end{bmatrix}
$$

となる。

対角ブロック以外が0であるのは、TDOA誤差とFDOA誤差が互いに独立であると仮定しているためである。

---

### 5.4 相関を持つ観測ノイズの生成

共分散行列 $R$ を持つ観測ノイズは、Cholesky分解を用いて生成できる。

共分散行列 $R$ を

$$
R
=
LL^T
$$

と分解する。

ここで、

- $L$：共分散行列 $R$ のCholesky因子である下三角行列
- $L^T$：行列 $L$ の転置行列
- Cholesky分解：正定値対称行列を下三角行列とその転置の積に分解する方法

である。

4次元の標準正規乱数ベクトルを

$$
\mathbf{n}
\sim
\mathcal{N}(\mathbf{0},I_4)
$$

とする。

ここで、

- $\mathbf{n}$：各成分が独立な標準正規分布に従う乱数ベクトル
- $\mathbf{0}$：すべての成分が0の平均ベクトル
- $I_4$：$4\times4$ 単位行列
- $\mathcal{N}(\mathbf{0},I_4)$：平均ベクトルが $\mathbf{0}$、共分散行列が $I_4$ の多変量正規分布

である。

この標準正規乱数ベクトルに $L$ を掛けて、

$$
\boxed{
\boldsymbol{\epsilon}
=
L\mathbf{n}
}
$$

とすると、

$$
\boldsymbol{\epsilon}
\sim
\mathcal{N}(\mathbf{0},R)
$$

となる。

すなわち、$\boldsymbol{\epsilon}$ は平均ベクトルが $\mathbf{0}$、共分散行列が $R$ の多変量正規分布に従う。

実際に共分散を計算すると、

$$
\begin{aligned}
\operatorname{Cov}(\boldsymbol{\epsilon})
&=
\operatorname{E}
\left[
\boldsymbol{\epsilon}
\boldsymbol{\epsilon}^T
\right]\\
&=
\operatorname{E}
\left[
L\mathbf{n}
\mathbf{n}^T L^T
\right]\\
&=
L
\operatorname{E}
\left[
\mathbf{n}\mathbf{n}^T
\right]
L^T\\
&=
LI_4L^T\\
&=
LL^T\\
&=
R
\end{aligned}
$$

となる。

ここで、

- $\operatorname{E}[\cdot]$：確率変数または確率ベクトルの期待値
- $\operatorname{Cov}(\cdot)$：共分散行列

を表す。

---

### 5.5 観測ノイズを含む観測モデル

理想的な観測ベクトルを $h(\mathbf{x})$、観測誤差ベクトルを $\boldsymbol{\epsilon}$ とすると、実際に得られる観測ベクトルは

$$
\boxed{
\mathbf{z}_{\rm obs}
=
h(\mathbf{x})
+
\boldsymbol{\epsilon}
}
$$

と表される。

ここで、

- $\mathbf{z}_{\rm obs}$：観測ノイズを含む実際の観測ベクトル
- $h(\mathbf{x})$：状態 $\mathbf{x}$ から計算される理想的な観測ベクトル
- $\boldsymbol{\epsilon}$：平均 $\mathbf{0}$、共分散行列 $R$ の観測誤差ベクトル

である。

観測誤差の確率分布は

$$
\boldsymbol{\epsilon}
\sim
\mathcal{N}(\mathbf{0},R)
$$

である。

したがって、実際の観測ベクトルは

$$
\mathbf{z}_{\rm obs}
\sim
\mathcal{N}
\left(
h(\mathbf{x}),
R
\right)
$$

と表すこともできる。

これは、観測ベクトルの平均値が $h(\mathbf{x})$、共分散行列が $R$ であることを意味する。