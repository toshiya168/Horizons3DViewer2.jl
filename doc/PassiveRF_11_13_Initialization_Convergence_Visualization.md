# 11. 初期値・12.収束判定・13.推定結果の可視化

# 11. 初期値の与え方

Gauss--Newton法およびLevenberg--Marquardt法は局所反復法であるため、初期値が収束性を大きく左右する。

初期状態を

$$
\mathbf{x}_0=
\begin{bmatrix}
x_0\\
y_0\\
v_{x0}\\
v_{y0}
\end{bmatrix}
$$

とする。

ここで

- $x_0,y_0$：初期位置
- $v_{x0},v_{y0}$：初期速度

である。

反復では

$$
\mathbf{x}_{k+1}
=
\mathbf{x}_k+\Delta\mathbf{x}
$$

を繰り返す。

初期値が真値に十分近ければ二次収束に近い高速収束を示すが、遠い場合は局所解への収束や発散が起こり得る。

本コードではシミュレーションで真値近傍を初期値として与え、アルゴリズム自体の性能を確認している。

---

# 12. 収束判定

各反復では更新量

$$
\Delta\mathbf{x}
$$

と残差

$$
\mathbf r=\mathbf z-h(\mathbf x)
$$

を評価する。

更新量のノルム

$$
\|\Delta\mathbf{x}\|
$$

が十分小さい場合、

$$
\|\Delta\mathbf{x}\|<\varepsilon_x
$$

で収束と判定する。

ここで

- $\varepsilon_x$：更新量の閾値

である。

また残差ノルム

$$
\|\mathbf r\|
$$

についても

$$
\|\mathbf r\|<\varepsilon_r
$$

なら十分観測を説明できたと判断できる。

さらに無限ループ防止のため

$$
k>N_{\max}
$$

となれば反復を終了する。

ここで

- $k$：反復回数
- $N_{\max}$：最大反復回数

である。

---

# 13. 推定結果の可視化

本サンプルでは位置平面と速度平面を同時表示する。

## 13.1 Position Plane

位置平面では

- 真値
- 推定位置
- TDOA双曲線
- Joint推定

を表示する。

TDOAのみでは

$$
\Delta t_{21},\;
\Delta t_{31}
$$

から得られる2本の双曲線の交点が位置推定となる。

---

## 13.2 Velocity Plane

速度平面では

- 真値速度
- 推定速度
- FDOA拘束線
- Joint推定速度

を表示する。

位置が既知なら

$$
A\mathbf v=b
$$

を解くことで速度が得られる。

Joint推定では位置も同時更新するため、

$$
\mathbf x=
[x,y,v_x,v_y]^T
$$

全体を一括推定する。

---

## 13.3 誤差の評価

位置誤差

$$
e_p
=
\left\|
\hat{\mathbf p}
-
\mathbf p_{\rm true}
\right\|
$$

速度誤差

$$
e_v
=
\left\|
\hat{\mathbf v}
-
\mathbf v_{\rm true}
\right\|
$$

を計算すると推定精度を評価できる。

ここで

- $\hat{\mathbf p}$：推定位置
- $\mathbf p_{\rm true}$：真値位置
- $\hat{\mathbf v}$：推定速度
- $\mathbf v_{\rm true}$：真値速度

である。

---

## 13.4 可視化の意義

位置平面ではTDOAが位置を拘束する様子、速度平面ではFDOAが速度を拘束する様子を直感的に確認できる。

さらにJoint推定では位置と速度が互いに影響し合いながら反復的に改善される過程を観察できるため、TDOA単独・FDOA単独・Joint推定の違いを理解しやすい。

# まとめ

本コードでは

1. 初期値を設定
2. LM反復で状態推定
3. 収束判定
4. Position PlaneとVelocity Planeで結果を可視化

という流れでPassive RFによる位置・速度推定を実現している。
