# TDOA/FDOA同時非線形推定において Joint Nonlinear Estimate と FDOA-only at Estimated Position が一致する理由

## 概要

未知状態

$$
\mathbf{x}=
\begin{bmatrix}
x\\y\\v_x\\v_y
\end{bmatrix}
$$

観測

$$
\mathbf{z}=
\begin{bmatrix}
\Delta t_{21}\\
\Delta t_{31}\\
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

観測数4、未知数4である。

## TDOA

$$
r_i=\|\mathbf p-\mathbf s_i\|
$$

$$
\Delta t_{21}=\frac{r_2-r_1}{c},
\qquad
\Delta t_{31}=\frac{r_3-r_1}{c}
$$

TDOAは位置のみの関数である。

## FDOA

LOSベクトル

$$
\mathbf u_i=
\frac{\mathbf p-\mathbf s_i}
{\|\mathbf p-\mathbf s_i\|}
$$

距離変化率

$$
\dot r_i=\mathbf u_i^T\mathbf v
$$

FDOA

$$
\Delta f_{21}
=
-\frac{f_0}{c}
(\mathbf u_2-\mathbf u_1)^T\mathbf v
$$

$$
\Delta f_{31}
=
-\frac{f_0}{c}
(\mathbf u_3-\mathbf u_1)^T\mathbf v
$$

位置が既知なら

$$
A\mathbf v=\mathbf b
$$

であり

$$
A=
\begin{bmatrix}
(\mathbf u_2-\mathbf u_1)^T\\
(\mathbf u_3-\mathbf u_1)^T
\end{bmatrix}
$$

$$
\mathbf b=
-\frac{c}{f_0}
\begin{bmatrix}
\Delta f_{21}\\
\Delta f_{31}
\end{bmatrix}
$$

より

$$
\mathbf v=A^{-1}\mathbf b
$$

となる。

## Joint Nonlinear Estimation

目的関数

$$
J=
\frac12
(\mathbf z-h(\mathbf x))^T
R^{-1}
(\mathbf z-h(\mathbf x))
$$

Gauss-Newton更新

$$
\delta=
(H^TR^{-1}H+\lambda I)^{-1}
H^TR^{-1}
(\mathbf z-h(\mathbf x))
$$

$$
\mathbf x_{k+1}=\mathbf x_k+\delta
$$

## 一致する理由

収束後の位置を

$$
\hat{\mathbf p}
$$

とすると、Joint推定の速度は

$$
A(\hat{\mathbf p})\mathbf v=\mathbf b
$$

を満たす。

一方、FDOA-onlyも同じ

$$
A(\hat{\mathbf p})\mathbf v=\mathbf b
$$

を直接解いている。

したがって

$$
\boxed{
\hat{\mathbf v}_{Joint}
=
\hat{\mathbf v}_{FDOA}
}
$$

となる。

## なぜか

現在は

- 観測4
- 未知数4

さらに

- TDOA→位置
- FDOA→速度

という構造なので

1. TDOAで位置決定
2. FDOAで速度決定

という逐次解法とJoint推定が同じ解へ収束する。

## 一致しなくなる条件

- 4局以上
- 複数時刻観測
- 加速度状態追加
- 観測数＞未知数

このときFDOAも位置推定へ寄与し、

$$
\hat{\mathbf v}_{Joint}
\neq
\hat{\mathbf v}_{FDOA}
$$

となる。

## 結論

現在の一致は偶然ではなく、3局・単一時刻・未知数4・観測4という問題設定の数学的帰結である。
