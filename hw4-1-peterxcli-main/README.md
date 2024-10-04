[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/SvA1ARLV)

# 多處理機平行程式設計 作業4-1說明


## 題目:讓專業的來



### 題目敘述

財經專家龍哥認為地球不是圓的，應該是橢圓的。所以不應該用方型矩陣來掃氣壓圖，而是要用一個 $D_1 × D_2$ 的矩陣$K$ 來掃氣壓圖，其中 $D_1$ 與 $D_2$ 均為奇數。令 $A_t[0..m − 1, 0..n − 1]$代表在 $t$ 時間時的氣壓圖矩陣，其中 $A_0$ 為題目給定的初始矩陣。下一個時間點的矩陣經由以下公式產生。



$A_{t+1} [i,j]$ $=$  $($ $\frac{1}{D_1 * D_2}$ $\sum_{\Delta_i= - \frac{D_1-1}{2}}^\frac{D_1-1}{2}$ $\space$  $\sum_{\Delta_j= - \frac{D_2-1}{2}}^\frac{D_2-1}{2}$ $\space$ $K [\frac{D_1-1}{2} + \Delta_i , \frac{D_2-1}{2} + \Delta_j]  A_t[i+\Delta_i,j+ \Delta_j]$       $)$ 






其中 $0$ $\leq$ $i$ $\leq$ $m − 1$$,$ $0$ $\leq$ $j$ $\leq$ $n − 1$。由於地球是圓的，所以當 $[i$ $+$ $\Delta$$i$ $,$ $j$ $+$ $\Delta$$i$$]$ 超過$[0..m − 1, 0..n − 1]$ 時需轉換成對應的座標，比如 $[−1, −1]$ 轉換成 $[m − 1, n − 1]、[−2, −2]$ 轉換成 $[m − 2, n − 2]、[m, n]$ 轉換成 $[0, 0]$ 以此類推。請使用  $pthread$ 來實作並輸出 $A_t$ 的結果。


### 輸入輸出說明

第一行只有一個數字 $t$，代表在 $t$ 時間的氣壓圖矩陣，第二行中的兩個數字  $n$ 和 $m$ 為氣壓圖矩陣的長寬。接下來$n$ $*$ $m$個數字是$A_0$裡的數值 $($$Row-Major$$)$，輸入完 $A_0$ 後下一行即是矩陣 $K$ 的長($D1$)和寬($D2$)以及矩陣$K$的內容。

例如：
![image](https://hackmd.io/_uploads/SkUtSVQNT.png)



資料範圍：
* $1 \leq n \leq 1000$
* $1 \leq m \leq 1000$
* $1 \leq t < 71$
* $1 \leq d_1 \leq 10$
* $1 \leq d_2 \leq 10$

輸出的 $A_t$ 內容請以$Row-Major$的方式印出，如以下格式
![image](https://hackmd.io/_uploads/ryJz8Vm4T.png)



<font color="#f00">注意:1.輸入的矩陣皆為int 
       $\space$$\space$$\space$$\space$$\space$$\space$$\space$ 2.輸出的數字後面都有空格</font>


### 繳交格式

在Github上傳一個程式碼檔案<font color="#f00">以及對應的Makefile檔案</font>，程式碼檔名為 學號_hw4_1 
例如：p12345678_hw4_1.c p12345678_hw4_1.cpp都可
Makefile就叫Makefile

