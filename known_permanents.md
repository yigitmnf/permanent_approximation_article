# Known Permanents Database

This file tracks matrices in `Matrices/` whose permanent is known or is
treated as ground truth for the current experiments.

Matrices with permanent `0` are intentionally omitted.
Families such as `Small`, `Medium`, `Large`, and `ErdosRenyi` are not listed
here because their permanents are currently unknown.

## TinyOriginal
Entries in this section come from `/home/kaya/Rasmussen/CleanBackup/PermanentResults/small_permanent`.

### Matrices/TinyOriginal/GD95_c.mtx
Permanent: 2295328928067623.5
Rows: 62
Cols: 62
Stored entries: 287
MatrixMarket: matrix coordinate pattern general
Matrix type: Directed graph from GD95
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/bcspwr01.mtx
Permanent: 376417000
Rows: 39
Cols: 39
Stored entries: 85
MatrixMarket: matrix coordinate pattern symmetric
Matrix type: Power network matrix
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/bcspwr02.mtx
Permanent: 17339123388
Rows: 49
Cols: 49
Stored entries: 108
MatrixMarket: matrix coordinate pattern symmetric
Matrix type: Power network matrix
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/chesapeake.mtx
Permanent: 13163539661940
Rows: 39
Cols: 39
Stored entries: 170
MatrixMarket: matrix coordinate pattern symmetric
Matrix type: Chesapeake Bay water quality model
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/curtis54.mtx
Permanent: 8.7765646985476352e+18
Rows: 54
Cols: 54
Stored entries: 291
MatrixMarket: matrix coordinate pattern general
Matrix type: Directed graph from Curtis
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/dwt_59.mtx
Permanent: 4.7271563194740265e+18
Rows: 59
Cols: 59
Stored entries: 163
MatrixMarket: matrix coordinate pattern symmetric
Matrix type: Discrete wavelet transform
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/ibm32.mtx
Permanent: 2398815
Rows: 32
Cols: 32
Stored entries: 126
MatrixMarket: matrix coordinate pattern general
Matrix type: IBM 32 processor matrix
Source: UF Sparse Matrix Collection

### Matrices/TinyOriginal/mycielskian6.mtx
Permanent: 1.3643529163476296e+18
Rows: 47
Cols: 47
Stored entries: 236
MatrixMarket: matrix coordinate pattern symmetric
Matrix type: Mycielski graph
Source: SuiteSparse Matrix Collection

### Matrices/TinyOriginal/will57.mtx
Permanent: 1069780304738061874
Rows: 57
Cols: 57
Stored entries: 281
MatrixMarket: matrix coordinate pattern general
Matrix type: Directed graph from Will
Source: UF Sparse Matrix Collection

## Banded
These permanents are exact and come from the cyclic band generator.

### Matrices/Banded/band_n1000_k2.mtx
Permanent: 104207051596066368588844988664977126343363937944778168928506449998330371001117189608703069071702510716775225113639692682101996584825295374597399059969286006359987999499408111127654711264867600511329765633181014532807761141236158080142135283307507082787222116035728466124162285581596160956440650882753501166613724691331758939774731632865784793498168836163183519941805840
Rows: 1000
Cols: 1000
Stored entries: 5000
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=1000, k=2
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

### Matrices/Banded/band_n1000_k3.mtx
Permanent: 939830720065424591189711712281170943756574697542470183834299424848606176195917799779209576525831926700110924764180435250349913171797902295565844402409699152690960135622409680140217230078035924558632742069739490265311941408959268399076856539328300948463239805458128001739368870470267968910234260605016974092727873478128937815339302960695953205126526330449409269706750997796049069204500819079340392207771862575251295242855131553598100329229410666193584948646734560579470734793642364297873
Rows: 1000
Cols: 1000
Stored entries: 7000
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=1000, k=3
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

### Matrices/Banded/band_n100_k2.mtx
Permanent: 6335628666511209600131664794318061520
Rows: 100
Cols: 100
Stored entries: 500
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=100, k=2
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

### Matrices/Banded/band_n100_k3.mtx
Permanent: 3956716713002412344334078800237720059846631839553
Rows: 100
Cols: 100
Stored entries: 700
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=100, k=3
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

### Matrices/Banded/band_n500_k2.mtx
Permanent: 10208185519281395220388226047624585105027636253377118204878889988543828549809452433197913114314902277194749720719101161797971718544244322459109728401823035594215331276860374037858093520
Rows: 500
Cols: 500
Stored entries: 2500
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=500, k=2
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

### Matrices/Banded/band_n500_k3.mtx
Permanent: 969448668092037037685333903236214416904649984997892773515940352867780953274679099492259344073156706923895428798699047467296025415730690875479539267682487004746428766839231313622010480094339770623945431101788499899869401271639590767586401364553
Rows: 500
Cols: 500
Stored entries: 3500
MatrixMarket: matrix coordinate integer general
Matrix type: cyclic band
Parameters: n=500, k=3
Construction: A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.

## BlockDiagonal
These permanents are exact because the matrices are block diagonal and each block permanent is computed exactly.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p1.mtx
Permanent: 3.8347441156680496167e+493
Rows: 1000
Cols: 1000
Stored entries: 6623
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=1000, block_size_counts={5:30, 6:32, 7:28, 8:24, 9:30}, reduction_ratio=0.1, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p1_weighted.mtx
Permanent: 1.4453882240130743876e+665
Rows: 1000
Cols: 1000
Stored entries: 6606
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=1000, block_size_counts={5:30, 6:32, 7:28, 8:24, 9:30}, reduction_ratio=0.1, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p3.mtx
Permanent: 2.620685437791253769e+390
Rows: 1000
Cols: 1000
Stored entries: 5303
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=1000, block_size_counts={5:30, 6:25, 7:36, 8:29, 9:24}, reduction_ratio=0.3, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p3_weighted.mtx
Permanent: 4.4471732910792458985e+571
Rows: 1000
Cols: 1000
Stored entries: 5397
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=1000, block_size_counts={5:30, 6:25, 7:36, 8:29, 9:24}, reduction_ratio=0.3, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p5.mtx
Permanent: 4.510363770646548555e+288
Rows: 1000
Cols: 1000
Stored entries: 4253
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=1000, block_size_counts={5:25, 6:27, 7:26, 8:36, 9:27}, reduction_ratio=0.5, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p5_weighted.mtx
Permanent: 2.7822789674401960788e+452
Rows: 1000
Cols: 1000
Stored entries: 4219
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=1000, block_size_counts={5:25, 6:27, 7:26, 8:36, 9:27}, reduction_ratio=0.5, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p7.mtx
Permanent: 4.0695442127306778373e+135
Rows: 1000
Cols: 1000
Stored entries: 2879
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=1000, block_size_counts={5:27, 6:36, 7:29, 8:31, 9:22}, reduction_ratio=0.7, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p7_weighted.mtx
Permanent: 9.712600710347080689e+302
Rows: 1000
Cols: 1000
Stored entries: 2889
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=1000, block_size_counts={5:27, 6:36, 7:29, 8:31, 9:22}, reduction_ratio=0.7, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p9.mtx
Permanent: 3.00578991243264e+16
Rows: 1000
Cols: 1000
Stored entries: 1623
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=1000, block_size_counts={5:22, 6:32, 7:24, 8:28, 9:34}, reduction_ratio=0.9, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n1000_p0p9_weighted.mtx
Permanent: 5.5356925392252322994e+176
Rows: 1000
Cols: 1000
Stored entries: 1633
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=1000, block_size_counts={5:22, 6:32, 7:24, 8:28, 9:34}, reduction_ratio=0.9, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p1.mtx
Permanent: 7.520305585124328796e+49
Rows: 100
Cols: 100
Stored entries: 674
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=100, block_size_counts={5:3, 6:2, 7:3, 8:2, 9:4}, reduction_ratio=0.1, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p1_weighted.mtx
Permanent: 6.156604212152369465e+66
Rows: 100
Cols: 100
Stored entries: 671
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=100, block_size_counts={5:3, 6:2, 7:3, 8:2, 9:4}, reduction_ratio=0.1, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p3.mtx
Permanent: 1.2131133049527156363e+42
Rows: 100
Cols: 100
Stored entries: 570
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=100, block_size_counts={5:5, 7:1, 8:4, 9:4}, reduction_ratio=0.3, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p3_weighted.mtx
Permanent: 1.03874766486162076154e+61
Rows: 100
Cols: 100
Stored entries: 580
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=100, block_size_counts={5:5, 7:1, 8:4, 9:4}, reduction_ratio=0.3, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p5.mtx
Permanent: 9.904970604918571463e+26
Rows: 100
Cols: 100
Stored entries: 416
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=100, block_size_counts={5:2, 6:3, 7:3, 8:3, 9:3}, reduction_ratio=0.5, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p5_weighted.mtx
Permanent: 1.4106026604282056956e+46
Rows: 100
Cols: 100
Stored entries: 429
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=100, block_size_counts={5:2, 6:3, 7:3, 8:3, 9:3}, reduction_ratio=0.5, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p7.mtx
Permanent: 42177191215104.0
Rows: 100
Cols: 100
Stored entries: 280
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=100, block_size_counts={5:3, 6:4, 7:4, 8:3, 9:1}, reduction_ratio=0.7, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p7_weighted.mtx
Permanent: 3.3658721583122436366e+30
Rows: 100
Cols: 100
Stored entries: 270
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=100, block_size_counts={5:3, 6:4, 7:4, 8:3, 9:1}, reduction_ratio=0.7, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p9.mtx
Permanent: 16.0
Rows: 100
Cols: 100
Stored entries: 155
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=100, block_size_counts={5:2, 6:4, 7:2, 8:2, 9:4}, reduction_ratio=0.9, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n100_p0p9_weighted.mtx
Permanent: 1.932929483245224609e+18
Rows: 100
Cols: 100
Stored entries: 157
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=100, block_size_counts={5:2, 6:4, 7:2, 8:2, 9:4}, reduction_ratio=0.9, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p1.mtx
Permanent: 8.6687000217437080496e+251
Rows: 500
Cols: 500
Stored entries: 3401
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=500, block_size_counts={5:13, 6:11, 7:17, 8:11, 9:18}, reduction_ratio=0.1, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p1_weighted.mtx
Permanent: 2.2232228375492189087e+339
Rows: 500
Cols: 500
Stored entries: 3416
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=500, block_size_counts={5:13, 6:11, 7:17, 8:11, 9:18}, reduction_ratio=0.1, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.1. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p3.mtx
Permanent: 2.9353293335637484008e+202
Rows: 500
Cols: 500
Stored entries: 2750
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=500, block_size_counts={5:8, 6:20, 7:11, 8:16, 9:15}, reduction_ratio=0.3, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p3_weighted.mtx
Permanent: 1.1481249628635831652e+290
Rows: 500
Cols: 500
Stored entries: 2739
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=500, block_size_counts={5:8, 6:20, 7:11, 8:16, 9:15}, reduction_ratio=0.3, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.3. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p5.mtx
Permanent: 4.9731284536925356646e+136
Rows: 500
Cols: 500
Stored entries: 2034
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=500, block_size_counts={5:16, 6:18, 7:17, 8:14, 9:9}, reduction_ratio=0.5, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p5_weighted.mtx
Permanent: 3.5247413005911163272e+226
Rows: 500
Cols: 500
Stored entries: 2045
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=500, block_size_counts={5:16, 6:18, 7:17, 8:14, 9:9}, reduction_ratio=0.5, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.5. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p7.mtx
Permanent: 2.503863712338368668e+66
Rows: 500
Cols: 500
Stored entries: 1450
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=500, block_size_counts={5:16, 6:12, 7:19, 8:10, 9:15}, reduction_ratio=0.7, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p7_weighted.mtx
Permanent: 3.771451799074192406e+148
Rows: 500
Cols: 500
Stored entries: 1427
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=500, block_size_counts={5:16, 6:12, 7:19, 8:10, 9:15}, reduction_ratio=0.7, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.7. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p9.mtx
Permanent: 1327104.0
Rows: 500
Cols: 500
Stored entries: 825
MatrixMarket: matrix coordinate integer general
Matrix type: reduced block diagonal
Parameters: n=500, block_size_counts={5:18, 6:9, 7:11, 8:18, 9:15}, reduction_ratio=0.9, weighted=False
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

### Matrices/BlockDiagonal/block_diagonal_n500_p0p9_weighted.mtx
Permanent: 4.7915297924512192044e+91
Rows: 500
Cols: 500
Stored entries: 827
MatrixMarket: matrix coordinate real general
Matrix type: reduced weighted block diagonal
Parameters: n=500, block_size_counts={5:18, 6:9, 7:11, 8:18, 9:15}, reduction_ratio=0.9, weighted=True, weight_range=[1, 2]
Construction: Diagonal entries are always kept. Off-diagonal entries inside each block are removed independently with probability 0.9. The block-diagonal matrix is then independently permuted in rows and columns. The exact permanent is the product of the exact block permanents.

## FromTinyToLarge
These permanents are products of the TinyOriginal values listed in this file.

### Matrices/FromTinyToLarge/generated_1000_971_5.mtx
Permanent: 52535155515535869128479212123383051841650589889186547732107923891577989184130328007677393209603820105114575796502940186821024180385875857890374665499777003253883805440578435340326367763354999393867752729766943057847210487637718364728661180416000000000000000000000000000
Rows: 971
Cols: 971
Stored entries: 5106
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x2, bcspwr01 x3, bcspwr02 x5, chesapeake x2, dwt_59 x3, ibm32 x1, mycielskian6 x3, will57 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_973_1.mtx
Permanent: 28458149807334034949063930033483443366626361231053559606942304859027995744527816789953889427824924517078066346583176460798985128396906364264791257860604759953847665195575897187200606258258371184010379653276016278258462723116212524142616036712417628458522256627529157915443200000000000000000000000000000
Rows: 973
Cols: 973
Stored entries: 5009
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x1, bcspwr02 x2, curtis54 x3, dwt_59 x6, ibm32 x1, mycielskian6 x2, will57 x3
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_973_3.mtx
Permanent: 1150990679808125021434609850931865040917703248535837769151406402874822097366494550200694332986869706069900848186692636410074621091669869949851342017808814253967837153120784630057972161346780580531352219515700794269861651655889398686111632642366207442771140608000000000000000000000000000000000
Rows: 973
Cols: 973
Stored entries: 5511
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x2, bcspwr01 x1, bcspwr02 x1, chesapeake x1, curtis54 x1, dwt_59 x5, ibm32 x4, mycielskian6 x4, will57 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_983_3.mtx
Permanent: 11369808642418517340007576541809315945535877080796919955327758426962075577247380576698671065851347298968370677431027296620059862385473515752800283001584550028347723489635486907462450458298876710310307459265380731696850712125075619840000000000000000000000000000000
Rows: 983
Cols: 983
Stored entries: 4609
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x3, bcspwr01 x4, bcspwr02 x3, chesapeake x1, curtis54 x3, dwt_59 x2, ibm32 x4, mycielskian6 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_989_4.mtx
Permanent: 343408552079758503307065712574141239882128370198060492440209272003428229789699898600815655095162911176417304174101534173806669444267362465409589121620471096607224603551032146381534124718377988645347788451034490218569274802620709057840244465649944887296000000000000000000000000000000
Rows: 989
Cols: 989
Stored entries: 4823
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x1, bcspwr01 x1, bcspwr02 x2, chesapeake x1, curtis54 x3, dwt_59 x4, ibm32 x6, mycielskian6 x1, will57 x2
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_993_4.mtx
Permanent: 114417182382229367609050061906517729913418137169043181159556270785284411487164705901878403669221175187707334697492552497533721513752343600247664681083572067895206241054049415937519110503000189987279635019248193287004098658038593345225906535021116514438742016000000000000000000000000000
Rows: 993
Cols: 993
Stored entries: 5253
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x2, bcspwr01 x3, bcspwr02 x4, chesapeake x3, curtis54 x1, dwt_59 x3, mycielskian6 x2, will57 x2
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_995_5.mtx
Permanent: 17688101594522588301615609552206937654409802636663365762981531386675855951687507750997684938154281630555109838538506690236584486764431195590893071566680738443787843123252112604446102566331067045938616475619751376230419386968620229035020399383599866682850847728134375840153600000000000000000000000
Rows: 995
Cols: 995
Stored entries: 5471
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x2, bcspwr01 x1, bcspwr02 x2, chesapeake x3, curtis54 x1, dwt_59 x3, ibm32 x2, mycielskian6 x2, will57 x4
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_1000_997_2.mtx
Permanent: 298028116286404189793746361021040509002482751697237120069536450817407143143535865645028497139579115442247761319824126672160018965590217403250594039605108707559258086463255166728963911795522906807502769726101174833648305751941205743974400000000000000000000000000000
Rows: 997
Cols: 997
Stored entries: 4494
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=1000
Composed from: GD95_c x3, bcspwr01 x5, bcspwr02 x2, curtis54 x1, dwt_59 x2, ibm32 x4, mycielskian6 x1, will57 x3
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_100_78_2.mtx
Permanent: 4954980108928468980000
Rows: 78
Cols: 78
Stored entries: 471
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=100
Composed from: bcspwr01 x1, chesapeake x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_100_88_3.mtx
Permanent: 6526740808340796000
Rows: 88
Cols: 88
Stored entries: 298
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=100
Composed from: bcspwr01 x1, bcspwr02 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_100_89_4.mtx
Permanent: 2566205041710233894279310
Rows: 89
Cols: 89
Stored entries: 407
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=100
Composed from: ibm32 x1, will57 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_100_96_1.mtx
Permanent: 14082095470981737181124762875560
Rows: 96
Cols: 96
Stored entries: 621
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=100
Composed from: chesapeake x1, will57 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_100_96_5.mtx
Permanent: 23656683561329191935721084800
Rows: 96
Cols: 96
Stored entries: 639
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=100
Composed from: bcspwr02 x1, mycielskian6 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_500_477_2.mtx
Permanent: 1502118503426628338462432220888431462721358613841406442730702369452443578462395842499248574983843363853103245521415014104170496000000000000
Rows: 477
Cols: 477
Stored entries: 2732
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=500
Composed from: GD95_c x2, bcspwr01 x1, bcspwr02 x1, chesapeake x2, curtis54 x2, ibm32 x1, mycielskian6 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_500_488_3.mtx
Permanent: 4763798086385994561665529334489176654717919575657376209188346413277341023139803970161405850727705636395034350586887666435961533347584000000000000
Rows: 488
Cols: 488
Stored entries: 2824
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=500
Composed from: GD95_c x1, bcspwr01 x1, bcspwr02 x1, chesapeake x1, dwt_59 x1, ibm32 x1, mycielskian6 x2, will57 x2
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_500_494_1.mtx
Permanent: 4711080986298265229322328984113891364912264883398371056099543054535359742464121412898054073587332970485078581665333170953115332760320000000000
Rows: 494
Cols: 494
Stored entries: 2544
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=500
Composed from: GD95_c x2, curtis54 x1, dwt_59 x1, ibm32 x3, mycielskian6 x1, will57 x2
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_500_497_5.mtx
Permanent: 1228478605989529987817799564536965179457327070088203227915353470555461483479581030485793614148861274672614566726897000646072598528000000000000000000000
Rows: 497
Cols: 497
Stored entries: 2919
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=500
Composed from: bcspwr01 x4, curtis54 x1, ibm32 x1, mycielskian6 x3, will57 x2
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

### Matrices/FromTinyToLarge/generated_500_499_4.mtx
Permanent: 41079282149586638101339366304921806597238843668471406094910144153459853827093555985102653004587250235322748940678340416315325600211756301680640000000
Rows: 499
Cols: 499
Stored entries: 2823
MatrixMarket: matrix coordinate integer general
Matrix type: block diagonal from TinyOriginal
Parameters: target_n=500
Composed from: GD95_c x4, chesapeake x1, curtis54 x2, mycielskian6 x1, will57 x1
Construction: Block diagonal composition of TinyOriginal matrices with independent row and column permutations.

