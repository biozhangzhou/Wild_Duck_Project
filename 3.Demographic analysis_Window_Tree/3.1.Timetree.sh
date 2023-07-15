#25-raml-5node.tree
25 1
((Ana-016,((((((PZ,HL),(LW,((ZW,(((CAU,SX),Ma),SB)),LC))),(Ana-005,Ana-003))'>.383587015703<10.65254126063',Ana-004),(((Ana-001,Ana-014)Ana-013),Ana-002))'>3.659030006337<24.14496428975',(((Ana-012,Ana-011),Ana-010),(Ana-009,((Ana-006,Ana-007),Ana-008))))'>6.492021576084<27.6370395394')'>23.14185735494<53.1465027444',Ana-015)'>30.76816305273<58.87333147957';
#mcmctree.ctl
###################################
         seed = -1
       seqfile = 2ALL-except-WZ-Filter-CAU-Chr-25-Bird.4d.fa.phy
      treefile = 25-raml-5node.tree
       outfile = 25sample

         ndata = 1
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
      * RootAge = '<10'  * constraint on root age, used if no fossil for root.

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5   * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1   * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 1047 1   * gamma prior for rate for genes   ### 1/替换率
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1  * times, rates, mixing, paras, RateParas

         print = 1
        burnin = 20000000
      sampfreq = 1000
       nsample = 50000

*** Note: Make your window wider (100 columns) when running this program.
####################################
#run "mcmctree mcmctree.ctl" with usedata=3
mv out.BV in.BV
#rerun "mcmctree mcmctree.ctl" with usedata=2
