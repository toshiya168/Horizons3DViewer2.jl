# 付録A コード全体の処理フロー

```text
真値の位置・速度設定
        │
        ▼
TDOA・FDOA観測生成
        │
        ▼
観測ノイズ付加
        │
        ▼
初期値設定
        │
        ▼
Levenberg–Marquardt反復
 ├─ 観測予測 h(x)
 ├─ 残差 r=z-h(x)
 ├─ 数値ヤコビアン H
 ├─ LM更新量 Δx
 └─ 状態更新
        │
        ▼
収束判定
        │
        ▼
位置・速度推定結果
        │
        ▼
Position Plane / Velocity Plane 可視化
```

---

# 付録B 主要関数一覧と役割

| 関数 | 役割 |
|------|------|
| `measurement_model(...)` | TDOA・FDOA観測値を計算する観測関数 |
| `numerical_measurement_jacobian(...)` | 中央差分により数値ヤコビアンを計算 |
| `levenberg_marquardt(...)` | LM法による非線形最小二乗推定 |
| `generate_measurement(...)` | 真値からノイズ付き観測を生成 |
| `add_noise(...)` | 共分散行列に基づき観測ノイズを生成 |
| `position_plane_plot(...)` | TDOA双曲線と推定位置を描画 |
| `velocity_plane_plot(...)` | FDOA拘束と推定速度を描画 |

これらの関数は以下のような依存関係となる。

```text
generate_measurement
        │
        ▼
measurement_model
        │
        ▼
numerical_measurement_jacobian
        │
        ▼
levenberg_marquardt
        │
        ▼
可視化関数
```

---

# 付録C 記号一覧（Notation）

| 記号 | 意味 |
|------|------|
| $\mathbf{x}$ | 状態ベクトル |
| $\mathbf{x}_0$ | 初期状態 |
| $\mathbf{p}$ | 位置ベクトル |
| $\mathbf{v}$ | 速度ベクトル |
| $\mathbf{z}$ | 観測ベクトル |
| $h(\mathbf{x})$ | 観測関数 |
| $\mathbf{r}$ | 残差ベクトル |
| $H$ | 観測ヤコビアン |
| $R$ | 観測ノイズ共分散 |
| $L$ | Cholesky分解の下三角行列 |
| $\sigma_t$ | TOA標準偏差 |
| $\sigma_f$ | 周波数標準偏差 |
| $r_i$ | 局$i$までの距離 |
| $\mathbf{u}_i$ | 局$i$方向の単位ベクトル |
| $c$ | 光速 |
| $f_0$ | 送信周波数 |
| $\Delta t_{ij}$ | TDOA |
| $\Delta f_{ij}$ | FDOA |
| $\Delta\mathbf{x}$ | 状態更新量 |
| $\lambda$ | LMダンピング係数 |
| $\Phi$ | 状態遷移（CR3BP等） |
| $I$ | 単位行列（単位写像） |
| $N_{\max}$ | 最大反復回数 |
| $\varepsilon_x$ | 更新量収束閾値 |
| $\varepsilon_r$ | 残差収束閾値 |

以上で、本資料で使用した主要な数式・アルゴリズム・記号を一通り網羅した。
