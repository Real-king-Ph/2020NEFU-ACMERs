# Math
## Number Theory

### 扩展欧几里得算法

__贝祖定理__  
&emsp;&emsp;若要使形如 $ax + by = c$ 的方程有解，那么必须要满足的条件是 $c = k\cdot\gcd(a,b)(k\in Z)$  

证明如下:  

首先不妨令 &emsp;&emsp;&emsp;&emsp;&emsp;$d = \gcd(a,b)$  
那么我们便得到 &emsp;&emsp;&emsp;$a'dx + b'dy = dc'$  
即 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;$a'x+b'y=c'$  
由于此时 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;$\gcd(a',b') = 1$  

1. 当$c' = 1$时  
所以我们得到的方程为$a'x + b'y = 1$  
由于当$a' = b'$时$a'=b'=1$  
不妨设 $a'>b'$  
由带余数除法我们可以得到，以下的方程  
$a'=b'\cdot q+r,0\leq r<b'$  
$b'=r\cdot q_1+r_1,0\leq r_1<r$  
$r=r_1\cdot q_2+r_2,0\leq r_2<r_1$  
$r_1=r_2\cdot q_3+r_3,0\leq r_3<r_2$  
$r_2=r_3\cdot q_4+r_4,0\leq r_4<r_3$  
……  
$r_{n-1}=r_n\cdot q_{n+1}+0$
>由于余数是单调递减的所以说我们的余数总会有等于 0 的时候。  
此时我们将上面化简的式子在一步一步的带回去就组成了我们原式的解。

2.  $c'\geq 2$ 时，只需在第一种情况下乘上系数即可。

那么写成程序方便写的形式

由 $ax+by=\gcd(a,b)$  
得 $a'x' + b'y' = \gcd(a,b)$  ，其中$a'= b,b'= a \mod b$  
由于 $a\mod b = a - \lfloor\frac{a}{b}\rfloor\times b$  
那么 $bx'+ (a- \lfloor\frac a b\rfloor \times b)y'=c$  
化简 $ay'+b(x'-\lfloor\frac a b \rfloor\times y') = c$  
那么 $x = y',b=x'-\lfloor\frac a b \rfloor\times y'$

```cpp
ll ExGCD(ll a, ll b,ll &x, ll &y){
   if(!b) {x = 1 , y = 0 ; return a ; }
   ll k = ExGCD(b,a % b, x, y);
   ll _x = x, _y = y;
   x = _y , y = _x - a / b * _y;
   return k;
}
```