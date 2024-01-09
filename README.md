
# PON-mt-tRNA
`PON-mt-tRNA` is a posterior probability-based method for classification of mitochondrial tRNA variations developed by Abhishek Niroula and Mauno Vihinen (PON-mt-tRNA: a multifactorial probability-based method for classification of mitochondrial tRNA variations Nucleic Acids Res. (2016) 44(5): 2020-2027.)[paper link](http://nar.oxfordjournals.org/content/early/2016/02/02/nar.gkw046.abstract)


# R shiny server for PON-mt-tRNA
We make a R shiny server for `PON-mt-tRNA`. 
You can run it locally (on your own computer) or public (make it available on the network). 


# Requirements 
- `Linux` with `docker` (If you want to make it available on the network).
- `R (4.2.3)` with `shiny(1.8.0)`, `shinydashboard(0.7.2)`, `DT(0.28)`, `digest(0.6.33)`.
- `Python (3.12.0)` with `pandas(2.1.4)`.
The versions we used are in brackets. Please note that the versions do not necessarily have to be the same. However, if you encounter an error, please check the version.


# Run PON-mt-tRNA locally
Clone this repository via the commands:
```  
git clone https://github.com/zhanghaoyang0/pon_mt_trna.git
cd pon_mt_trna
Rscript shiny.r
```

Once the above has completed, you will get a local website for PON-mt-tRNA: 

You will see: 
<video controls>
  <source src="https://github.com/zhanghaoyang0/pon_mt_trna/blob/master/www/video.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>


# Run PON-mt-tRNA public (on the network)
If you want to deposit PON-mt-tRNA in a port (e.g., 8503 in the example) of your own (public-access) server (e.g., xxx.com): 
``` 
docker pull zhanghaoyang0/rshiny
docker run -itd -p 8503:3838 -v path_of_pon_mt_trna:/srv/shiny-server/pon_mt_trna zhanghaoyang0/rshiny
``` 
Then, you can visit PON-mt-tRNA in a port of your server (e.g., xxx.com:8503).


# Feedback and comments
Add an issue or send email to haoyang.zhang@med.lu.se