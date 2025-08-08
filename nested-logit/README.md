このフォルダに置かれているプログラムは、近視眼的なNested Logitモデルに従い行動する旅行者の観光周遊行動モデルの推定と予測を行います。

旅行者は観光地を訪問する度に、以下のNested Logitモデルに従い選択を行います。
上位ネストでは「出国vs周遊継続」の選択を行います．周遊継続の選択確率は以下の式で与えられます．

```math
P_\mathrm{continue} = \frac { \exp \left( \gamma V_\mathrm{continue} \right) } { \exp \left( \gamma V_\mathrm{continue} \right) + \exp \left( - \gamma \beta \bar{c}_{\mathrm{current},\mathrm{departure}} \right)  }
```
```math
V_\mathrm{continue} = \log \sum_{r \in K_{\mathrm{unvisited}}} \exp \left( \alpha_r - \beta \bar{c}_{\mathrm{cur},r} \right)
```

```math
P_k = \frac { \exp \left( \alpha_k - \beta \bar{c}_{\mathrm{cur},k} \right) } { \sum_{r \in K_{\mathrm{unvisited}}} \exp \left( \alpha_r - \beta \bar{c}_{\mathrm{cur},r} \right)  }
```
