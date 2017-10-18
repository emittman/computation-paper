library(drake)

dir.create("arrayOutput")
# set.seed(170915123)
# sample(1e9, 4)
#[1] 836728172 355907037 978695416 548741190

my_plan <- plan(
  inits1 = initialize_chain(seed = 836728172, K = 2^11, n.iter = 100000, methodPi = "stickBreaking"),
  inits2 = initialize_chain(seed = 355907037, K = 2^11, n.iter = 100000, methodPi = "stickBreaking"),
  inits3 = initialize_chain(seed = 978695416, K = 2^11, n.iter = 100000, methodPi = "stickBreaking"),
  inits4 = initialize_chain(seed = 548741190, K = 2^11, n.iter = 100000, methodPi = "stickBreaking"),
  chain1_sb = sample_bnp_model(inits1),
  chain2_sb = sample_bnp_model(inits2),
  chain3_sb = sample_bnp_model(inits3),
  chain4_sb = sample_bnp_model(inits4),
  strings_in_dots = "literals"
)

check(my_plan)
# plot_graph(my_plan)
save(my_plan, file = "data/targets.RData")