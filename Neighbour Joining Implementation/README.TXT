This program generates a  tree using neighbor-joining method.

Requirements: Python 2.7 (http://www.python.org/) , NumPy (http://numpy.scipy.org/)
Usage:
  python nj.py <input file>

#<input file>: (space-separated distance matrix)

Example:  
  python nj.py test_data.txt

R code for plot can be found in 'R codeForPlot.txt'
Requirements: R (version 3.0.2 and above)

Installation instructions:
http://rpy.sourceforge.net/rpy2/doc-dev/html/overview.html
https://cran.r-project.org/
https://launchpad.net/~marutter/+archive/ubuntu/rrutter

https://www.youtube.com/watch?v=Nxl7HDUyw0I  (Video tutorial to get latest version: 3.3.2)

#Calculating RF Distance between Biopython NJ tree and my NJ implementation:
Usage:
  python dendropyWithImplementedNJ.py

Observations:

My NJ Tree (Newick):
(lcl|AP004411.1_gene_4:227.4128418,((lcl|AP008989.1_gene_3:101.46875,lcl|AP008990.1_gene_3:87.53125):35.1015625,lcl|AP008988.1_gene_3:177.3984375):12.0871582,((((((((((((lcl|AM489718.1_gene_16:3.911764706,lcl|AB182305.1_gene_3:3.088235294):21.74038462,lcl|AM489716.1_gene_16:25.25961538):3.958333333,(lcl|AM919429.1_gene_16:26.16071429,lcl|NC_010121.1_gene_4:25.83928571):4.166666667):2.176136364,lcl|NC_012323.1_gene_3:32.51136364):27.734375,(lcl|NC_015097.1_gene_16:38.68333333,lcl|FR751399.1_gene_16:36.31666667):27.328125):1.081597222,lcl|NC_007396.1_gene_3:70.16840278):33.01367188,lcl|FR751401.1_gene_16:97.87695312):13.453125,lcl|AP004412.1_gene_3:98.0234375):13.42415365,lcl|NC_015120.1_gene_16:117.1266276):14.58764648,lcl|NC_002616.1_gene_3:175.3908691):5.958170573,lcl|AP008991.1_gene_3:181.1980794):8.954589844,lcl|AP002929.1_gene_3:151.2080078):1.041992188);

Biopython NJ Tree:
((lcl|AP002929.1_gene_3:0.09916,(lcl|AP008991.1_gene_3:0.11755,(lcl|NC_002616.1_gene_3:0.11297,(lcl|NC_015120.1_gene_16:0.07565,(lcl|AP004412.1_gene_3:0.06929,(lcl|FR751401.1_gene_16:0.06286,(lcl|NC_007396.1_gene_3:0.04507,((lcl|NC_012323.1_gene_3:0.02088,((lcl|NC_010121.1_gene_4:0.01660,lcl|AM919429.1_gene_16:0.01680)Inner4:0.00268,(lcl|AM489716.1_gene_16:0.01622,(lcl|AB182305.1_gene_3:0.00198,lcl|AM489718.1_gene_16:0.00251)Inner1:0.01396)Inner5:0.00254)Inner6:0.00140)Inner7:0.01781,(lcl|FR751399.1_gene_16:0.02332,lcl|NC_015097.1_gene_16:0.02484)Inner3:0.01755)Inner8:0.00069)Inner9:0.02120)Inner10:0.00937)Inner11:0.00819)Inner12:0.00905)Inner14:0.00409)Inner15:0.00483)Inner16:0.00280,((lcl|AP008990.1_gene_3:0.05646,lcl|AP008989.1_gene_3:0.06878)Inner2:0.02216,lcl|AP008988.1_gene_3:0.11432)Inner13:0.00680,lcl|AP004411.1_gene_4:0.14895)Inner17:0.00000;

Distance between Biopython NJ tree and my NJ Implementation tree:
0

    
