# Passive RF測位の基礎（1～3）

## 1. Passive RF測位とは

Passive RF（Passive Radio Frequency）測位は、自ら電波を送信せず、対象物から放射される電波を複数の受信局で受信し、対象物の**位置**と**速度**を推定する技術である。

ここでは、対象物と受信局が二次元平面上に存在する場合を考える。

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

- $x,y$：対象物の位置座標
- $v_x,v_y$：対象物の速度成分

である。

対象物の位置ベクトルを

$$
\mathbf{p}
=
\begin{bmatrix}
x\\
y
\end{bmatrix}
$$

とし、対象物の速度ベクトルを

$$
\mathbf{v}
=
\begin{bmatrix}
v_x\\
v_y
\end{bmatrix}
$$

とする。

また、受信局 $i$ の位置ベクトルを

$$
\mathbf{s}_i
=
\begin{bmatrix}
s_{i,x}\\
s_{i,y}
\end{bmatrix}
$$

とする。

ここで、

- $i$：受信局を識別する番号
- $s_{i,x},s_{i,y}$：受信局 $i$ の位置座標

である。

対象物と受信局 $i$ の間の距離を $r_i$ とすると、

$$
r_i
=
\|\mathbf{p}-\mathbf{s}_i\|
$$

である。

ここで、$\|\cdot\|$ はベクトルのユークリッドノルムを表す。したがって、

$$
r_i
=
\sqrt{
(x-s_{i,x})^2
+
(y-s_{i,y})^2
}
$$

となる。

Passive RF測位では、主に次の観測量を利用する。

- TDOA（Time Difference of Arrival：到来時刻差）
- FDOA（Frequency Difference of Arrival：到来周波数差）

---

## 2. TDOAの原理

受信局 $i$ における電波の受信時刻を $t_i$ とすると、

$$
t_i
=
t_0+\frac{r_i}{c}
$$

である。

ここで、

- $t_i$：受信局 $i$ における受信時刻
- $t_0$：対象物が電波を送信した時刻
- $r_i$：対象物と受信局 $i$ の間の距離
- $c$：光速

である。

送信時刻 $t_0$ は一般に未知であるため、複数の受信局間で受信時刻の差を考える。

受信局1を基準局とした、受信局2と受信局1の受信時刻差を

$$
\Delta t_{21}
=
t_2-t_1
$$

と定義する。

ここで、記号 $\Delta$ は二つの量の差を表し、添字 $21$ は「受信局2の値から受信局1の値を引く」ことを表す。

受信時刻の式を代入すると、

$$
\Delta t_{21}
=
\left(
t_0+\frac{r_2}{c}
\right)
-
\left(
t_0+\frac{r_1}{c}
\right)
$$

となるため、

$$
\Delta t_{21}
=
\frac{r_2-r_1}{c}
$$

が得られる。

同様に、受信局3と受信局1の受信時刻差を

$$
\Delta t_{31}
=
t_3-t_1
$$

と定義すると、

$$
\Delta t_{31}
=
\frac{r_3-r_1}{c}
$$

となる。

TDOAに光速を掛けることで、距離差に変換できる。

受信局2と受信局1に対する距離差を

$$
\Delta r_{21}
=
r_2-r_1
$$

と定義すると、

$$
\Delta r_{21}
=
c\Delta t_{21}
=
r_2-r_1
$$

となる。

したがって、対象物の位置 $\mathbf{p}$ は

$$
\|\mathbf{p}-\mathbf{s}_2\|
-
\|\mathbf{p}-\mathbf{s}_1\|
=
\Delta r_{21}
$$

を満たす。

この式を満たす点の集合は、受信局1と受信局2を焦点とする**双曲線**である。

同様に、

$$
\Delta r_{31}
=
r_3-r_1
=
c\Delta t_{31}
$$

と定義すると、

$$
\|\mathbf{p}-\mathbf{s}_3\|
-
\|\mathbf{p}-\mathbf{s}_1\|
=
\Delta r_{31}
$$

が得られる。

したがって、2本の双曲線

$$
\|\mathbf{p}-\mathbf{s}_2\|
-
\|\mathbf{p}-\mathbf{s}_1\|
=
\Delta r_{21}
$$

および

$$
\|\mathbf{p}-\mathbf{s}_3\|
-
\|\mathbf{p}-\mathbf{s}_1\|
=
\Delta r_{31}
$$

の交点が、対象物の位置候補となる。


---

## 3. FDOAの原理

対象物が移動している場合、対象物と受信局の間の距離が時間とともに変化するため、ドップラー効果によって受信周波数が変化する。

受信局 $i$ で観測される受信周波数を $f_i$ とすると、低速近似のもとで

$$
f_i
=
f_0
\left(
1-\frac{\dot r_i}{c}
\right)
$$

と表される。

ここで、

- $f_i$：受信局 $i$ で観測される受信周波数
- $f_0$：対象物が送信する電波の送信周波数
- $\dot r_i$：対象物と受信局 $i$ の間の距離速度
- $c$：光速

である。

記号 $\dot r_i$ は、距離 $r_i$ の時間微分を表す。

$$
\dot r_i
=
\frac{dr_i}{dt}
$$

ここで、$t$ は時間である。

$\dot r_i>0$ の場合は対象物と受信局の距離が増加しており、$\dot r_i<0$ の場合は距離が減少している。

受信局 $i$ から対象物へ向かう視線方向の単位ベクトルを $\mathbf{u}_i$ とすると、

$$
\mathbf{u}_i
=
\frac{\mathbf{p}-\mathbf{s}_i}
{\|\mathbf{p}-\mathbf{s}_i\|}
$$

で定義される。

ここで、

- $\mathbf{u}_i$：受信局 $i$ から対象物へ向かう視線単位ベクトル
- $\mathbf{p}$：対象物の位置ベクトル
- $\mathbf{s}_i$：受信局 $i$ の位置ベクトル

である。

視線単位ベクトルは二次元の場合、

$$
\mathbf{u}_i
=
\begin{bmatrix}
u_{i,x}\\
u_{i,y}
\end{bmatrix}
$$

と表される。

受信局が静止している場合、距離速度は

$$
\dot r_i
=
\mathbf{u}_i^T\mathbf{v}
$$

となる。

ここで、

- $\mathbf{v}$：対象物の速度ベクトル
- $\mathbf{u}_i^T$：視線単位ベクトル $\mathbf{u}_i$ の転置
- 上付きの $T$：行列またはベクトルの転置

を表す。

成分表示すると、

$$
\dot r_i
=
u_{i,x}v_x
+
u_{i,y}v_y
$$

である。

したがって、$\dot r_i$ は対象物の速度 $\mathbf{v}$ を視線方向 $\mathbf{u}_i$ に射影した速度成分である。

受信局2と受信局1の周波数差を

$$
\Delta f_{21}
=
f_2-f_1
$$

と定義する。

ここで、$\Delta f_{21}$ は受信局2で観測された周波数から、受信局1で観測された周波数を引いた値である。

受信周波数の式を代入すると、

$$
\Delta f_{21}
=
f_0
\left(
1-\frac{\dot r_2}{c}
\right)
-
f_0
\left(
1-\frac{\dot r_1}{c}
\right)
$$

となるため、

$$
\Delta f_{21}
=
-\frac{f_0}{c}
(\dot r_2-\dot r_1)
$$

が得られる。

さらに、

$$
\dot r_1
=
\mathbf{u}_1^T\mathbf{v},
\qquad
\dot r_2
=
\mathbf{u}_2^T\mathbf{v}
$$

を代入すると、

$$
\boxed{
\Delta f_{21}
=
-\frac{f_0}{c}
(\mathbf{u}_2-\mathbf{u}_1)^T\mathbf{v}
}
$$

となる。

同様に、受信局3と受信局1の周波数差を

$$
\Delta f_{31}
=
f_3-f_1
$$

と定義すると、

$$
\boxed{
\Delta f_{31}
=
-\frac{f_0}{c}
(\mathbf{u}_3-\mathbf{u}_1)^T\mathbf{v}
}
$$

となる。

対象物の位置 $\mathbf{p}$ が既知であれば、各視線単位ベクトル $\mathbf{u}_i$ を計算できる。そのため、2つのFDOA観測式は、速度ベクトル $\mathbf{v}$ に関する連立一次方程式として

$$
A\mathbf{v}
=
\mathbf{b}
$$

と表せる。

ここで、

- $A$：視線単位ベクトルの差から構成される観測行列
- $\mathbf{b}$：FDOA観測値から構成されるベクトル

である。

具体的には、

$$
A
=
\begin{bmatrix}
(\mathbf{u}_2-\mathbf{u}_1)^T\\
(\mathbf{u}_3-\mathbf{u}_1)^T
\end{bmatrix}
$$

および

$$
\mathbf{b}
=
-\frac{c}{f_0}
\begin{bmatrix}
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

である。

二次元の場合、

$$
A
=
\begin{bmatrix}
u_{2,x}-u_{1,x} & u_{2,y}-u_{1,y}\\
u_{3,x}-u_{1,x} & u_{3,y}-u_{1,y}
\end{bmatrix}
$$

となる。

行列 $A$ が正則、すなわち

$$
\det(A)\neq 0
$$

であれば、速度ベクトルは

$$
\boxed{
\mathbf{v}
=
A^{-1}\mathbf{b}
}
$$

として求められる。

ここで、

- $\det(A)$：行列 $A$ の行列式
- $A^{-1}$：行列 $A$ の逆行列

である。

## 4 傾きの一致
---

受信局1と受信局2の観測によるTDOA双曲線を表す関数を

$$
g_{21}(\mathbf{p})
=
\|\mathbf{p}-\mathbf{s}_2\|
-
\|\mathbf{p}-\mathbf{s}_1\|
-
\Delta r_{21}
$$

と定義する。

TDOA双曲線上では、

$$
g_{21}(\mathbf{p})=0
$$

が成立する。

位置ベクトル $\mathbf{p}$ に関する距離 $r_i$ の勾配は、

$$
\nabla_{\mathbf{p}} r_i
=
\frac{\mathbf{p}-\mathbf{s}_i}
{\|\mathbf{p}-\mathbf{s}_i\|}
=
\mathbf{u}_i
$$

である。

ここで、$\nabla_{\mathbf{p}}$ は位置ベクトル $\mathbf{p}$ に関する勾配を表す。

したがって、TDOA双曲線を表す関数 $g_{21}$ の勾配は、

$$
\nabla_{\mathbf{p}}g_{21}
=
\mathbf{u}_2-\mathbf{u}_1
$$

となる。

勾配ベクトルは、その点における双曲線の接線に垂直な法線ベクトルである。したがって、TDOA双曲線上のある点における微小変位を

$$
d\mathbf{p}
=
\begin{bmatrix}
dx\\
dy
\end{bmatrix}
$$

とすると、接線方向では

$$
(\mathbf{u}_2-\mathbf{u}_1)^T d\mathbf{p}
=
0
$$

が成立する。

成分表示すると、

$$
(u_{2,x}-u_{1,x})dx
+
(u_{2,y}-u_{1,y})dy
=
0
$$

である。

したがって、$u_{2,y}-u_{1,y}\neq 0$ の場合、受信局1と受信局2によるTDOA双曲線の接線の傾きは、

$$
\boxed{
\frac{dy}{dx}
=
-
\frac{u_{2,x}-u_{1,x}}
{u_{2,y}-u_{1,y}}
}
$$

となる。

同様に、受信局1と受信局3によるTDOA双曲線の接線の傾きは、

$$
\boxed{
\frac{dy}{dx}
=
-
\frac{u_{3,x}-u_{1,x}}
{u_{3,y}-u_{1,y}}
}
$$

となる。

後述するように、対象物の位置 $\mathbf{p}$ を固定した場合、受信局1と受信局2によるFDOA制約は速度平面上で

$$
(\mathbf{u}_2-\mathbf{u}_1)^T\mathbf{v}
=
-\frac{c}{f_0}\Delta f_{21}
$$

という直線になる。

このFDOA直線の傾きは、

$$
\frac{dv_y}{dv_x}
=
-
\frac{u_{2,x}-u_{1,x}}
{u_{2,y}-u_{1,y}}
$$

である。

したがって、同じ対象位置 $\mathbf{p}$ で視線単位ベクトルを評価する限り、

$$
\boxed{
\text{TDOA双曲線の接線の傾き}
=
\text{FDOA直線の傾き}
}
$$

が厳密に成立する。

同様に、受信局1と受信局3の組についても、

$$
\boxed{
\text{TDOA}_{31}\text{双曲線の接線の傾き}
=
\text{FDOA}_{31}\text{直線の傾き}
}
$$

が成立する。

したがって、対象位置における2本のTDOA双曲線の接線交差角と、同じ位置で構成した2本のFDOA直線の交差角も、厳密に一致する。

この一致は遠方近似によるものではない。TDOA双曲線の位置に関する勾配と、FDOA観測式の速度に関する係数が、いずれも

$$
\mathbf{u}_j-\mathbf{u}_1
$$

という同じ視線単位ベクトル差で表されることに由来する。