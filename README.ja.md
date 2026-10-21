<h1 align='center'>sharpInterFoam</h1>

[![CI](https://github.com/sugurushiratori/sharpInterFoam/actions/workflows/ci.yml/badge.svg)](https://github.com/sugurushiratori/sharpInterFoam/actions/workflows/ci.yml)

THINC/SW 法と S-CLSVOF 法を組み合わせた、非構造格子対応の OpenFOAM 用二相流ソルバーです。
次の論文のコードです:
[*Int. J. Microgravity Sci. Appl.*, 38(3) 380301 (2021)](https://doi.org/10.15011/jasma.38.380301)

[English README](README.md)

## ソルバー

| ソルバー | 内容 |
|---|---|
| `sharpInterFoam` | 非圧縮・等温の二相流（interFoam に、体積率の THINC/SW 移流と、界面の法線・曲率を求める S-CLSVOF のレベルセット関数を組み込んだもの） |
| `sharpInterThermoFoam` | `sharpInterFoam` にエネルギー式と表面張力の温度依存性（マランゴニ対流・熱毛管流）を加えたもの |

主な機能

* 非構造格子上の THINC/SW による界面移流（向き付き流束・陽的更新・サブサイクル）
* S-CLSVOF：体積率からレベルセット関数を構成（atanh による初期値と、narrow band に限定した Godunov 型の再初期化）
* 表面張力は CSF 法または density-scaled CSF 法（圧力勾配と釣り合う離散化）
* 混合密度は体積率、またはレベルセットの平滑化 Heaviside 関数から評価（質量と運動量の輸送が整合）
* 表面張力の温度依存性 `sigma(T) = sigma + sigmaTemp*(T - TRef)` を法線方向と接線方向（マランゴニ力）の両方に反映
* 毛管時間刻み制約、接触角、移動格子、並列計算（領域分割によらず同じ結果）

## 計算例

| 上昇気泡（[Hysing et al. 2009](validation/risingBubble)） | Zalesak 円盤（[例題](tutorials/ZalesakDisc)） |
|:---:|:---:|
| <img src="docs/images/risingBubble.gif" height="320"> | <img src="docs/images/zalesakDisc.gif" height="320"> |

<p align="center"><img src="docs/images/thermocapillaryCavity.gif" width="560"><br>
熱毛管キャビティ：温度・流線・界面（<a href="tutorials/thermocapillaryCavity">例題</a>）</p>

## 対応する OpenFOAM

OpenFOAM（ESI / OpenCFD）**v2412** と **v2512**。

## ビルド

```sh
source <OpenFOAM のインストール先>/etc/bashrc
./Allwmake          # ライブラリ libTHINC, libSCLSVOF と 2 つのソルバー
./Allwclean         # ビルド結果の削除
```

ライブラリとソルバーは `$FOAM_USER_LIBBIN` と `$FOAM_USER_APPBIN` に置かれます。

## 使い方

例題は [`tutorials`](tutorials) にあります。

| ケース | ソルバー | 内容 |
|---|---|---|
| `ZalesakDisc/structured`, `ZalesakDisc/unstructured` | sharpInterFoam | 切り欠き円盤の剛体回転（論文 4.1 節） |
| `dropletDeformation/structured` | sharpInterFoam | 静止液滴とスプリアス流れ（論文 4.2 節） |
| `thermocapillaryCavity` | sharpInterThermoFoam | キャビティ内の熱毛管流（論文 4.3 節, 400 x 160） |

各ケースは `./Allrun`（並列計算は `./AllrunPar`、あるケースのみ）で実行し、`./Allclean` で初期状態に戻します。

### パラメータ

`system/fvSolution` の体積率の辞書:

```
"alpha.water.*"
{
    nAlphaSubCycles 2;

    SCLSVOF                          // レベルセット関数
    {
        correctPsi      yes;         // yes: S-CLSVOF; no: interFoam と同じく alpha から曲率を計算
        widthFactor     1.5;         // 界面の半幅 epsilon = widthFactor × 格子幅
        densityScaled   yes;         // density-scaled CSF モデル（密度比が大きい場合）
        densityFunctionHeaviside yes;// 混合密度を Heaviside 関数（yes）または alpha（no）から計算
        initializeAtanh yes;         // 初期値 psi0 = epsilon/pi*atanh(2 alpha - 1)（yes）または線形（no）
        denomDeltaTau   10;          // 再初期化の擬似時間刻み = 格子幅 / denomDeltaTau
        factorNumLoop   1;           // 再初期化の反復回数 = widthFactor*denomDeltaTau*factorNumLoop
        tolMagGradPsi   1e-8;        // 再初期化の残差の相対許容値
    }

    THINC
    {
        enable          yes;         // yes: THINC/SW; no: MULES（interFoam）
        epsAlpha        1e-12;       // alpha < epsAlpha または alpha > 1 - epsAlpha では風上流束
        epsBeta         1e-2;        // THINC の急峻さの下限。これより小さい範囲は風上流束と線形に混合
    }
}
```

`constant/transportProperties`:

```
phases (water air);              // alpha = 1 が water
water { transportModel Newtonian; nu 1e-6; rho 1000; }
air   { transportModel Newtonian; nu 1.5e-5; rho 1; }
sigma   0.07;

// sharpInterThermoFoam: 両相とも transportModel NewtonianThermo とし、
// 温度拡散率 kappa [m^2/s] と、必要なら Cp [J/(kg K)] を与える
// （Cp がない場合は rho*Cp = 1、熱伝導率 = kappa）
sigmaTemp   -1.5e-4;             // d sigma / d T
TRef        293;                 // sigma は温度 TRef での表面張力（既定値 0）
```

`system/controlDict`:

```
adjustTimeStep  yes;
maxCo           0.25;
maxAlphaCo      0.25;
maxCapillaryNo  1;               // dt <= maxCapillaryNo * sqrt((rho1 + rho2) h^3 / (4 pi sigma))
```

### 推奨設定

* 時間スキームは `Euler`（Crank-Nicolson では質量と運動量の輸送の整合が崩れます）。
* `adjustTimeStep yes` とし、`maxCapillaryNo` は 1 以下。
* `interpolationSchemes { default linear; }`。`pointLinear` ではスプリアス流れが大きくなります。
* `nOuterCorrectors` を大きくして `residualControl` で打ち切ると、レベルセット再初期化の計算量が外部反復の回数倍になります。まずは 1〜数回で試してください。

## テストと検証

* [`tests/regression`](tests/regression): CI で実行する小規模な回帰テスト（`tests/regression/run.sh`）。
* [`validation`](validation): 評価スクリプト付きのベンチマーク。結果は [`docs/ja/validation.md`](docs/ja/validation.md)。
* [`tests/transposeCheck`](tests/transposeCheck): 一部の OpenFOAM ビルドの不具合の検査（下記）。

## 既知の問題

* 一部の OpenFOAM ビルドでは、一時テンソル場の転置 `T()` が誤った値を返します
  （conda-forge 版 OpenFOAM v2412, gcc 13 で確認）。この場合、3 次元計算と
  x-z 面の 2 次元計算では、どの OpenFOAM ソルバーでも層流の粘性応力が誤ります。
  お使いの環境は `tests/transposeCheck/Allrun` で確認できます。

## 論文版からの変更

[CHANGELOG.md](CHANGELOG.md) を参照してください。論文で使用した版にはタグ `v1.0-ijmsa2021` を付けています。

## 引用

```bibtex
@article{Shiratori_IJMSA_2021,
  title   = {Efficient Implementation of Two-Phase Flow Solver Based on {THINC/SW} and {S-CLSVOF} on Unstructured Meshes},
  author  = {Shiratori, Suguru and Usui, Takuro and Koyama, Shiho and Ozawa, Shumpei and Nagano, Hideaki and Shimano, Kenjiro},
  journal = {International Journal of Microgravity Science and Application},
  volume  = {38},
  number  = {3},
  pages   = {380301},
  year    = {2021},
  doi     = {10.15011/jasma.38.380301}
}
```

## ライセンス

GNU General Public License v3（[COPYING.txt](COPYING.txt)）。sharpInterFoam は OpenFOAM の派生物であり、
OpenCFD Ltd. および OpenFOAM Foundation の承認を受けたものではありません。
