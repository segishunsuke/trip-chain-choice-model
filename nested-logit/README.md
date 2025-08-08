このフォルダに置かれているプログラムは、以下のNested Logitモデルに従い行動する旅行者の観光周遊行動モデルの推定と予測を行います。

旅行者は観光地を訪問する度に、次に訪問する
```math
P_k = \frac { \exp \left( \alpha_k - \beta \bar{c}_{\mathrm{cur},k} \right) } { \sum_{r \in K_{\mathrm{unvisited}}} \exp \left( \alpha_r - \beta \bar{c}_{\mathrm{cur},r} \right)  }
```
