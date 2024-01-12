
# R shiny server for PON-mt-tRNA
`PON-BTK` is a method for classifying variations in the kinase domain of Bruton tyrosine kinase (BTK) related to X-linked agammaglobulinemia (XLA) to disease-causing and harmful.  

The method is using the semi-supervised classification method ['upclass'](https://cran.r-project.org/web/packages/upclass/index.html). It uses an expectation maximization (EM) algorithm to obtain maximum likelihood estimates of the model parameters and classifications for the unlabeled data.  

We developed a R shiny server for `PON-BTK`. You can visit the server at [here](http://lap676.srv.lu.se:8503/pon_btk/). The code and data of this server are public available. 

Reference: Jouni Väliaho, Imrul Faisal, Csaba Ortutay, C. I. Edvard Smith and Mauno Vihinen.
Characterization of all possible single nucleotide change-caused amino acid substitutions in the kinase domain of Bruton tyrosine kinase.
*Hum Mutat*. 2015. [paper link](https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22791).  


# Requirements 
- `Linux`.
- `R (4.2.3)` with `shiny(1.8.0)`, `shinydashboard(0.7.2)`, `DT(0.28)`.
The versions we used are in brackets. Please note that the versions do not necessarily have to be the same.


# Run the server locally
Clone this repository via the commands:
```  
git clone https://github.com/zhanghaoyang0/pon_btk.git
cd pon_btk
Rscript shiny.r
```
Then you will get a local website for PON-BTK.


# Feedback and comments
Add an issue or send an email to haoyang.zhang@med.lu.se