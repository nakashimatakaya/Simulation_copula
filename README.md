# Simulation study using Copula

## 研究の目的

RCTの実データからConditional Average Treatment Effect (CATE)の推定値を求め，

RCTのTable 1に記載された要約統計量（平均や標準偏差）から周辺分布を推定し，

copulaを用いたシミュレーションで得られたCATE推定値の妥当性を示す

## 研究の進め方

1.  保有するRCTの実データを用いてCATE値の推定をおこなうのに加えて，そのデータから各変数（BMI，採血結果など）の要約統計量（平均や分散）を算出．

2.  その後，その要約統計量から，各変数の周辺分布を仮定し，vine copulaを用いてn次元同時分布を推定．

3.  推定されたn次元同時分布から，各変数をsimulateし，simulateされたデータを用いてCATEを推定します．(CATEはCausal forestなどのモデルを用いて推定予定です．)

4.  最後に，実データに基づくCATEとシミュレートされたデータに基づくCATEを比較検討します．

## 参考にした論文，サイト

-   [Virtual Patient Simulation Using Copula Modeling](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/cpt.3099?campaign=wolearlyview)

    -   公開されていたgithub code：<https://github.com/vanhasseltlab/copula_vps>

-   [Introduction to R package copulaSim](https://cran.r-project.org/web/packages/copulaSim/vignettes/introduction.html)

-   SASだけど．．[Simulate multivariate correlated data by using PROC COPULA in SAS](https://blogs.sas.com/content/iml/2021/07/07/proc-copula-sas.html)
