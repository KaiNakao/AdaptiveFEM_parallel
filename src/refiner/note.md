refinement
- edge collapsing
- edge splitting
- face to edge [flip23]
- edge to face [flip32]
- edge swapping [flip44]
- node movement

境界近傍でメッシュ改善が発生すると決め打ちして考えてみるか

指定された要素から求まるShteiner点を追加してマーキング領域を制約付きDelaunay分割．


waiting list:
- マーキング
- 連なっている要素集合を拾ってくる
- シュタイナー点の計算
- Delaunay分割
- 曲面の復元

うわーたいへん！