library(drake)

# set.seed(170915123)
# sample(1e9, 4)
#[1] 836728172 355907037 978695416 548741190

my_plan <- plan(
  inits1 = initialize_chain(836728172, "stickBreaking"),
  inits2 = initialize_chain(355907037, "stickBreaking"),
  chain1_sb = sample_bnp_model(inits1),
  chain2_sb = sample_bnp_model(inits2),
  strings_in_dots = "literals"
)

check(my_plan)
build_graph(my_plan)

make(plan = my_plan, verbose = TRUE, jobs=2, packages = c("cudarpackage","dplyr"), )