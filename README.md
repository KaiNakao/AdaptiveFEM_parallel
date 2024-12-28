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
以下のジョブスクリプトが整備されている．
- job_gen_initial_mesh.sh
    - 地形データから初期メッシュを生成
    - WORKDIR, NPROCには，ワーキングディレクトリ名，順解析時のMPIプロセス数をそれぞれ設定する．
- job_forward_analysis.sh
    - 順解析を実行
    - WORKDIR, NPROCには，ワーキングディレクトリ名，順解析時のMPIプロセス数をそれぞれ設定する．
- job_refinement.sh
    - メッシュのrefinementを実行
    - WORKDIR, NPROCには，ワーキングディレクトリ名，順解析時のMPIプロセス数をそれぞれ設定する．
- job_gen_refined_mesh.sh
    - refinementの結果から順解析用のrefined meshを生成
    - WORKDIR_ORG, WORKDIR_NEW, NPROCには，refinementを行ったワーキングディレクトリ名，refined meshを格納するワーキングディレクトリ名，順解析時のMPIプロセス数をそれぞれ設定する．