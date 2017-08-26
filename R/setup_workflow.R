#' Constructs a dataframe identifying targets
#' Use drake to check dependency structure

library(drake)
dir.create("arrayOutput")

# set.seed(170915123)
# sample(1e9, 4)
#[1] 836728172 355907037 978695416 548741190

my_plan <- plan(
  inits1 = initialize_chain(836728172, "stickBreaking"),
  inits2 = initialize_chain(355907037, "stickBreaking"),
  inits3 = initialize_chain(978695416, "stickBreaking"),
  inits4 = initialize_chain(548741190, "stickBreaking"),
  chain1_sb = sample_bnp_model(inits1),
  chain2_sb = sample_bnp_model(inits2),
  chain3_sb = sample_bnp_model(inits3),
  chain4_sb = sample_bnp_model(inits4),
  strings_in_dots = "literals"
)

check(my_plan)
# plot_graph(my_plan)
save(my_plan, "targets.RData")