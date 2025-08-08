このフォルダのプログラムは、近視眼的Nested Logitモデルに基づく旅行者の観光周遊行動の推定と予測を行います。

旅行者は観光地を訪問する度に、二層のNested Logitモデルに従い選択を行います。
上位ネストでは「周遊終了」と「周遊継続」の選択を行います。周遊継続の選択確率は次式で与えられます。

```math
P_\mathrm{continue} = \frac { \exp \left( \gamma V_\mathrm{continue} \right) } { \exp \left( \gamma V_\mathrm{continue} \right) + \exp \left( - \gamma \beta \bar{c}_{\mathrm{current},\mathrm{departure}} \right)  }
```
```math
V_\mathrm{continue} = \log \sum_{l \in K_{\mathrm{unvisited}}} \exp \left( \alpha_l - \beta \bar{c}_{\mathrm{current},l} \right)
```
ただし、現在地がトリップチェイン始点の空港の場合には、上記の式は適用されず、必ず周遊継続が選択されます。以上の式において、$`\mathrm{current}`$は現在地、$`\mathrm{departure}`$はトリップチェイン終点の空港、$`K_{\mathrm{unvisited}}`$は未訪問の観光地の集合です。$`\gamma`$は$`0 < \gamma \le 1`$を満たすスケールパラメータです。

周遊継続が選択された場合は、下位ネストで次に訪問する観光地の選択を行います。未訪問の観光地$`k`$の選択確率は以下の式で与えられます．
```math
P_k = \frac { \exp \left( \alpha_k - \beta \bar{c}_{\mathrm{current},k} \right) } { \sum_{l \in K_{\mathrm{unvisited}}} \exp \left( \alpha_l - \beta \bar{c}_{\mathrm{current},l} \right)  }
```
