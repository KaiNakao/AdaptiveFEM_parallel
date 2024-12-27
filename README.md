## Adaptive FEM
ibisでの実行を想定．
前準備として，jivsmディレクトリで地形データをダウンロードしておく．
```sh
wget https://www.jishin.go.jp/main/chousa/12_choshuki/dat/nankai/lp2012nankai-e_str.zip 

unzip lp2012nankai-e_str.zip
```
workディレクトリでジョブを実行する．
```sh
qub **.sh
```
ジョブスクリプトのWORKDIR，NPROCには，
ワーキングディレクトリ名，順解析時のMPIプロセス数をそれぞれ設定する．
- job_gen_initial_mesh.sh
    - 地形データから初期メッシュを生成
- job_forward_analysis.sh
    - 順解析を実行
- job_refinement.sh
    - メッシュのrefinementを実行