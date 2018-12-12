# DMRG
Density Matrix Re-normalization Group and its variants

<br>This is on going java project for DMRG methods
<br>The project create under Netbeans IDEA
<br>I use common math library from apache :
<br>http://commons.apache.org/proper/commons-math/download_math.cgi
<br>make sure you download add the library in your project before you use the program

Initial ver:0.1
<br>Currently I have 2 runnable class, HeisenbergChain.java and DMRG.java
<br>Eigen vectors and eigen values calculate used Janczos and Davidson methods
<br>eigen value valid untill 16 particles, more than that result "run out memory"

Future work
1. Create dynamic local file to restore temporary data such alpha, beta, psy matrix to save memory space
2. The DMRG wave-function : Matrix product State (MPS) [available soon]
3. The density matrix renormalization group for ab initio quantum chemistry
