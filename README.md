# Implementation-of-NETS: Extremely Fast Outlier Detection from a Data Stream via Set-Based Processing

This is the implementation of the paper published in PVLDB 2019 [[Paper](http://www.vldb.org/pvldb/vol12/p1303-yoon.pdf)]

## 1. Overview
REMOVE AND REPLACE: This paper addresses the problem of efficiently detecting outliers from a data stream as old data points expire from and new data points enter the window incrementally. The proposed method is based on a newly discovered characteristic of a data stream, that the change in the locations of data points in the data space is typically very insignificant. This observation has led to the finding that the existing distance-based outlier detection algorithms perform excessive unnecessary computations that are repetitive and/or canceling out the effects. Thus, in this paper, we propose a novel set-based approach to detecting outliers, whereby data points at similar locations are grouped and the detection of outliers or inliers is handled at the group level. Specifically, a new algorithm NETS is proposed to achieve a remarkable performance improvement by realizing set-based early identification of outliers or inliers and taking advantage of the "net effect" between expired and new data points. Additionally, NETS is capable of achieving the same efficiency even for a high-dimensional data stream through two dimensional-level filtering.  Comprehensive experiments using six real-world data streams show 5 to 25 times faster processing time than state-of-the-art algorithms with comparable memory consumption. We assert that NETS opens a new possibility to real-time data stream outlier detection. Our goal for this implementation is to show with code adjustments and proper testing that NETS can perform even better than originally recorded.

## 2. Data Sets
| Name    | # data points  | # Dim    | Size    | Link           |
| :-----: | :------------: | :------: |:-------:|:--------------:|
| GAU     | 1M             | 1        |  7.74MB  |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/gaussian.txt) |
| STK     | 1.05M          | 1        |  7.57MB |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/stock.txt) |
| TAO     | 0.58M          | 3        |  10.7MB |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/tao.txt) |
<!--| HPC     | 1M             | 7        |  28.4MB  |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/household2.txt) |
| GAS     | 0.93M          | 10       |  70.7MB  |[link](http://archive.ics.uci.edu/ml/machine-learning-databases/00362/HT_Sensor_UCIsubmission.zip) |
| EM      | 1M             | 16       |  119MB  |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/ethylene.txt) |
| FC      | 1M             | 55       |  72.2MB  |[link](https://infolab.usc.edu/Luan/Outlier/Datasets/fc.data) |-->

## 3. Configuration
We are currently trying to implement the NETS algorithm in JAVA and run on **JDK 17.0.4.1.**
Code from utilities and non algorithm functions was brought over
Certain printouts, variables, and constructors were eliminated to improve time
Main algorithms were recoded

## 4. Proposed Parameters
```
--dataset: title of datasets (string, one of {GAU, STK, TAO, HPC, GAS, EM, FC})
--W: the size of a window (integer)
--S: the size of a slide (integer)
--R: the distance threshold (double)
--K: the number of neighbors threshold (integer)
--D: the number of full dimensions (integer)
--sD: the number of sub dimensions (integer)
--nW: the number of windows (integer)

## 5. Compilation and Running:
From the root directory of the project run: javac .\test\testBase.java
Sample usage after compilation: java test.testBase --dataset TAO --W 10000 --S 500 --R 1.9 --K 50 --D 3 --sD 3 --nW 1
Parameters for running can be changed accordingly
Testing is done on Florida State Universities Computer Science Departments servers
