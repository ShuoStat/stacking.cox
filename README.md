# Using Stacking to Integrate Clinical and Omics Data for Survival Prediction

In this study, we proposed a stacking-based approach for integrative predictions, introducing two strategies: component-wise base-learner stacking (Stk(CP)) and component-wise penalty factor base-learner stacking (Stk(PF)). The key difference between them lies in how base learners are fitted. These strategies were implemented for survival prediction and compared against competitors using multiple real-world datasets. The results demonstrated the superior accuracy and robustness of stacking-based methods. Furthermore, an extension of Stk(PF), termed Stk(SS), showcased the flexibility of stacking by enabling customized base learners for diverse datasets.



### Project Structure

```
.
├── codes/
├── data/
├── output/
├── results/
└── README.md
```

`codes` codes for reproducing.  
`data` datasets used for model building.  
`output` intermediate outputs.    
`results` final results.

### Files in codes

`cox.fun.R` functions for training models  
`drfs.R` dimension reduction based on feature screening   
`ipflasso` IPFlasso
`stk.glm_v4.R` and `stk.cox_v7` stacking method
`drpca.R` dimension reduction based on PCA  
`helper.R` helper function

### Reproducing the Results

Run *mainFunction.R* step by step
```
R /codes/mainFunction.R
```

### Contact

Maintainer: Shuo Wang  
GitHub: https://github.com/ShuoStat/  
Email: wangsures@foxmail.com  

### Custom License

BY-NC-ND 4.0 License © 2025 Shuo Wang

This work is provided for personal and educational use only.

You may not:
- Use this work in any form of academic publication, including but not limited to journal articles, conference papers, or theses.
- Use this work for any commercial purposes, including but not limited to products, services, or monetized platforms.
- Distribute modified or derivative versions of this work.

All rights are reserved. For inquiries, contact wangsures@foxmail.com
