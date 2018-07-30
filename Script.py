
# coding: utf-8

# In[1]:

from Generator import *
from FD_ADMM_Workers import *


# # Generate Problem instances and save  in data folder.

# In[ ]:


#Saves Instance_treek.pkl for k in [0,...,number_of_instances]
number_of_instances = 5
for i in range(number_of_instances):
    create_problem(graph = h,graph_type = 16, n_requests = 100000, n_split=  4, alpha_list = [1,2,4], n_cores = [1,2,4,16,32], name = i)


# # Computes Optimum values and save in  data folder

# In[ ]:

optims(number_of_instances)
#computes optimum for all and saves Optim_treek.pkl for all k


# # Runs on instances and save in data folder

# In[ ]:

#Choose the deadline in seconds; solves Instance_Treename, and saves as Tree_Resultname.pkl
run_one(name,deadline)

# Format: (Gps,Hst,It)

#It:
#dict with alpha = 1,2,4 as keys
#For each alpha, list of tuples (n_cores,N) where n_cores is the number of cores used and N is the number of achieved iterations

#Hst:
#dict with alpha = 1,2,4 as keys
#For each alpha, list of tuples (n_cores, W) where W is a dict in format {time1: optimality_gap1, timex: optimality_gapx, ...}

#Gps:
#dict with alpha = 1,2,4 as keys
#For each alpha, list of tuples (n_cores, G) where G is the achieved optimality gap.

