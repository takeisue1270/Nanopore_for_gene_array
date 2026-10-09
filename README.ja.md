# A Targeted Nanopore Sequencing Pipeline for Tandem Gene Array Engineering（タンデム遺伝子アレイ解析用の標的 Nanopore シークエンス解析パイプライン）

*[English README](README.md)*

## 概要

タンデム遺伝子アレイとは、同一の配列単位が染色体上で連続して並んだ構造です。
遺伝子重複を用いた手法で、この単位数を増減させた株を作製できます。

本リポジトリは、そのような株を標的 Nanopore シークエンスで解析するパイプラインです。
解析はバーコードごとに行います。アレイの単位数と構造上の特徴を評価します。
処理の内容は、リードのアライメント、標的配列と側方配列の抽出、カバレッジの算出、
および `.fsa` のリファレンスを用いた正規化です。

本体は `Alignment.sh` です。すべての手順はコマンドライン引数と環境変数で指定します。
実行条件を記録すれば、同じ処理を再現できます。

正規化（手順 5）は出芽酵母（*Saccharomyces cerevisiae*）のみに対応します。
染色体名が `chrI` 形式のリファレンスから各染色体の長さを読む実装のためです。
これ以外の手順は生物種に依存しません。

本パイプラインは研究用です。[ライセンス](LICENSE)に記載のとおり、無保証で提供します。
診断、治療、臨床上の判断には使用できません。遺伝子組換え生物を扱う実験は法令の規制を
受けます。日本ではカルタヘナ法が該当します。実施の可否は所属機関の安全委員会が判断します。

---

## 処理の流れ

バーコードごとに次の順で処理します。「任意」と付した手順は、環境変数で実行の
有無を指定します。

1. 最小リード長で絞り込みながら、バーコードごとに FASTQ を結合
2. **minimap2** でリファレンスにアライメントし、ソート済み BAM を作成
3. WIG のカバレッジトラックを作成（任意）
4. bedgraph のカバレッジトラックを作成（任意）
5. `.fsa` のリファレンスで bedgraph を正規化（任意、出芽酵母専用）
6. 標的のタンデムアレイ全体にまたがるリードを抽出（任意）
7. 標的にヒットしたリードのみを抽出（任意）
8. 標的の上流・下流の側方配列を抽出（任意）
9. 抽出した配列をアライメントし、YASS の SVG を作成（任意）
10. BLASTN による確認・分類（任意）

---

## 必要なもの

### 本体のツール

* minimap2
* samtools
* bedtools
* seqkit
* python 3.8 以降

### 任意のツール

* igvtools（WIG の作成に使用）
* yass（ドットプロットの作成に使用）
* BLAST+（BLASTN に使用）

### Python のパッケージ

```bash
pip install mappy biopython numpy
```

* `mappy` — minimap2 の Python バインディング。3 つの抽出スクリプトが使用します
* `biopython`・`numpy` — `Bedgraph_normalize.py` のみが使用します

---

## 各自で用意するリファレンス

次のファイルを `./NSA/` に置きます。いずれもリポジトリには含めていません。
内容が標的と使用するアセンブリによって異なるためです。また、この大きさの配列
ファイルは Git での管理に適しません。

| ファイル | 内容 | 用意のしかた |
| --- | --- | --- |
| `${REFBASE}.fa` | アライメントのリファレンス。改変したアレイを扱う場合は、ゲノム全体ではなくコンストラクトやアレイの 1 単位を指定するのが通例です。 | 配列編集ソフトから FASTA で書き出します。minimap2 の `.mmi` は初回に自動で作られます。 |
| `${REFBASE}.UP1000.fa` | 標的のアレイの直上流 1 kb。アレイ全体にまたがるリードの検出に使用します。 | 標的領域の配列から切り出します。名前は `REFBASE` と揃えます。 |
| `${REFBASE}.DOWN1000.fa` | 直下流 1 kb。用途は同じです。 | 同上。 |
| `ACT1.fa` | 単一コピーの対照領域。`RUN_BLAST=1` で検体ごとの分母として使用します。 | 使用する株の *ACT1* の配列です。単一コピーの領域であれば他でも同じ役割を果たします。 |
| `S288C_reference_sequence_R64-2-1_20150113_0CUP1RU_1rDNARU_phiX.fsa` | 出芽酵母 S288C R64-2-1 のリファレンスを改変したものです。CUP1 の繰り返し単位を 0、rDNA の繰り返し単位を 1 とし、phiX を加えています。`RUN_FSA_NORM=1` と側方配列のマッピングで使用します。 | [SGD の R64-2-1](https://www.yeastgenome.org/) を元に、対象の繰り返し単位について同じ改変を行うか、`FSA_REF` に各自のリファレンスを指定します。正規化を行うには染色体名を `chrI` 形式のままにしてください。 |

正規化用のリファレンスで繰り返し単位を減らすのは、意図した処理です。繰り返しを
リファレンスのコピー数のまま残した場合、改変したアレイのカバレッジ比は解釈でき
ません。どの単位をどの数に変えたかを記録してください。この記録がない正規化後の
値は解釈できません。

---

## 実行

```bash
./Alignment.sh PREFIX START END REFBASE [THREADS] [MINLEN]
```

実行例を次に示します。

```bash
./Alignment.sh T975 1 96 S288C 14 200
```

引数と環境変数の一覧は [README.md](README.md) にあります。出力はすべて実行時の
カレントディレクトリに作られます。同名のファイルは上書きされます。

---

## 引用

本パイプラインの引用には、次の Zenodo のアーカイブを用います。

* Takesue H. A Targeted Nanopore Sequencing Pipeline for Tandem Gene Array Engineering. Zenodo. doi:[10.5281/zenodo.18440611](https://doi.org/10.5281/zenodo.18440611)

パイプラインが呼び出すツールも併せて引用します。その一覧は [README.md](README.md) の
Citation にあります。`RUN_FSA_NORM=1` を使用した場合は、次の正規化スクリプトの
アーカイブも引用します。

* Okada S. poccopen/Bedgraph_norm_ratio. Zenodo. doi:[10.5281/zenodo.11515695](https://doi.org/10.5281/zenodo.11515695)

---

## 他者による著作物

`NSA/Bedgraph_normalize.py` は、Satoshi Okada 氏（[@poccopen](https://github.com/poccopen)）が
作成した **Bedgraph_norm_ratio** の一部です。本リポジトリでは、これを CC BY 4.0 の
条件で再配布しています。このファイルは本リポジトリの MIT ライセンスの対象では
ありません。条件の詳細は [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) にあります。

---

## ライセンス

本リポジトリのコードは MIT ライセンスです。条文は [LICENSE](LICENSE) にあります。
`NSA/Bedgraph_normalize.py` は例外で、前項のとおり CC BY 4.0 です。

Oxford Nanopore Technologies、ONT、MinION、GridION、PromethION は
Oxford Nanopore Technologies plc の商標です。本リポジトリで言及しているその他の
ツールの権利は、それぞれの著作者に帰属します。本リポジトリはこれらの開発元とは
無関係です。承認も支援も受けていません。名称は、読み込む対象と実行するツールを
示す目的で使用しています。
