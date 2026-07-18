

# 8. 数値ヤコビアン

Gauss--Newton法では観測関数のヤコビアン

$$
H=
\frac{\partial\mathbf h}{\partial\mathbf x}
$$

が必要となる。

ここで

- $H$：観測ヤコビアン行列
- $\mathbf h$：観測関数
- $\mathbf x$：状態ベクトル

本コードでは解析微分ではなく**中央差分**で近似している。

状態変数 $x_j$ を微小量 $\delta_j$ だけ変化させる。

$$
\mathbf x^{+}
=
\mathbf x+\delta_j\mathbf e_j
$$

$$
\mathbf x^{-}
=
\mathbf x-\delta_j\mathbf e_j
$$

ここで

- $\mathbf e_j$：第$j$成分のみ1の単位ベクトル
- $\delta_j$：微小変位

である。

すると偏微分は

$$
\boxed{
\frac{\partial\mathbf h}{\partial x_j}
\approx
\frac{
\mathbf h(\mathbf x^{+})
-
\mathbf h(\mathbf x^{-})
}
{2\delta_j}
}
$$

で近似できる。

これを全状態について計算すると

$$
H=
\begin{bmatrix}
\frac{\partial h_1}{\partial x}&
\frac{\partial h_1}{\partial y}&
\frac{\partial h_1}{\partial v_x}&
\frac{\partial h_1}{\partial v_y}\\
\frac{\partial h_2}{\partial x}&
\cdots&\cdots&\cdots\\
\vdots&&&\\
\frac{\partial h_4}{\partial x}&
\cdots&\cdots&\cdots
\end{bmatrix}
$$

となる。

コードでは

```julia
numerical_measurement_jacobian(...)
```

がこの計算を担当し、

```julia
(z_plus-z_minus)/(2*step)
```

によって中央差分を実装している。

中央差分は前進差分

$$
\frac{h(x+\delta)-h(x)}{\delta}
$$

よりも打切り誤差が小さく、

$$
O(\delta^2)
$$

の精度を持つため、本コードではこちらが採用されている。
