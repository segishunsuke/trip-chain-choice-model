このフォルダに置かれているプログラムは、以下のNested Logitモデルの推定と予測を行います。

```math
P_k = \frac { \exp \left( \alpha_k - \beta \bar{c}_{\mathrm{cur},k} \right) } { \sum_{r \in K_{\mathrm{unvisited}}} \exp \left( \alpha_r - \beta \bar{c}_{\mathrm{cur},r} \right)  }
```
