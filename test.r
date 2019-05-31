
library(e1071)
getwd()
install.packages("here")
source("./r/includes/libraries.r")

glimpse(mrgdf)
dandf <- mrgdf %>%
  select(ftype, RECID, wt, c05800, c05800.reform)
write_csv(dandf, "d:/tcdir/danrecs2.csv")

mrgdf %>%
  group_by(ftype) %>%
  summarise_at(vars(c05800, c05800.reform), ~sum(. * wt) / 1e6) %>%
  mutate(diff=.[[3]] - .[[2]],
         pdiff=diff / .[[2]] * 100)




# run amt on the full regular puf
fbase <- read_csv(paste0("d:/tcdir/", "puf-13-#-#-#.csv")) %>% mutate(ftype="base")
fref <- read_csv(paste0("d:/tcdir/", "puf-13-#-amtex_2013-#.csv")) %>% mutate(ftype="reform")
df <- bind_rows(fbase, fref)
glimpse(df)

df2 <- df %>% 
  pivot_longer(-c(ftype, RECID, FLPDYR, s006)) %>%
  pivot_wider(names_from = ftype) %>%
  mutate(diff=reform - base)

df2 %>%
  group_by(name) %>%
  summarise(n=n(), nup=sum(reform>base), ndown=sum(base>reform), nunchanged=sum(base==reform),
            nwtdret=sum(s006), nwtdup=sum((reform>base)*s006), nwtddown=sum((base>reform)*s006), nwtdunchanged=sum((base==reform)*s006),
            valchange=sum(diff * s006)) %>%
  kable(digits=0, format.args = list(big.mark=",")) %>% 
  kable_styling()

dfw <- df %>%
  pivot_wider(names_from = ftype, values_from = c(c00100, c09600, c62100, taxbc))
glimpse(dfw)

dfw %>% mutate(amtti_diff=c62100_reform - c62100_base,
               amtl_diff=c09600_reform - c09600_base,
               taxbc_diff=taxbc_reform - taxbc_base) %>%
  summarise(n_tax_increase=sum(taxbc_diff>0), n_tax_decrease=sum(taxbc_diff<0), n_nochange=sum(taxbc_diff==0))


x# AMT mystery
# first, get the baseline tax calc results for the relevant files
synth10 <- readRDS(paste0(globals$tc.dir, "synthpop10", "_rwprep.rds"))
names(synth10)

syn20 <-readRDS(paste0(globals$tc.dir, "synpuf20", "_rwprep.rds"))
names(syn20)

# get the stacked file on which reforms were run so that we have ftype and RECID for the reforms
test_stack <- readRDS(paste0(globals$tc.fir, "test_stack.rds"))

reform <- "amtex_2013"
refname <- paste0("synth10syn20_", reform, "_output.rds")
refdf <- readRDS(paste0(globals$tc.dir, refname))
names(refdf) %>% sort

# get the big merged file of the reform
refdf2 <- left_join(test_stack %>% dplyr::select(RECID, ftype),
                    refdf,
                    by="RECID") %>%
  mutate(reform_name=reform)
names(refdf2) %>% sort

# construct a file for analysis
pufbase <- synth10$tc.base %>%
  filter(ftype=="puf") %>%
  left_join(synth10$tc.output %>% 
              select(RECID, taxbc.base=taxbc,
                     c62100.base=c62100,
                     c09600.base=c09600), by=c("RECID"))

pufref <- refdf2 %>%
  filter(ftype=="puf") %>%
  select(RECID, taxbc.reform=taxbc, c62100.reform=c62100, c09600.reform=c09600)

pufcomp <- pufbase %>%
  left_join(pufref, by="RECID")

amtx <- 80800
amtxpostart <- 153900
amtxpor <- .25
(amtxpoend <- amtxpostart + amtx / amtxpor)

comp <- pufcomp %>%
  filter(MARS==2, c62100.base >= amtxpostart, c62100.base <= amtxpoend) %>%
  mutate(amtti.diff=c62100.reform - c62100.base,
         amtl.diff=c09600.reform - c09600.base,
         tax.diff=taxbc.reform - taxbc.base) %>%
  arrange(-abs(amtti.diff))

write_csv(comp, paste0(globals$tc.dir, "comp.csv"))

diffs <- comp %>%
  mutate(diff=taxbc.reform - taxbc.base) %>%
  filter(diff!=0)

diffs <- comp %>%
  mutate(diff=c62100.reform - c62100.base) %>%
  filter(diff!=0)


# amt variables
# Output Variable Name: c09600
# Description: Alternative Minimum Tax (AMT) liability
# 
# Output Variable Name: c62100
# Description: Alternative Minimum Tax (AMT) taxable income
# 
# tc Name: AMT_CG_rt1
# TB Name: Long term capital gain and qualified dividends (AMT) rate 1
# Description: Capital gain and qualified dividends (stacked on top of regular income) below threshold 1 are taxed at this rate.























library(profvis)
library(dplyr)

profvis({
  diamonds <- read.csv("optimization/data/diamonds.csv")
  
  diamonds_by_cut <- diamonds %>%
    group_by(cut) %>%
    summarise_if(is.numeric, mean)
  
  write.csv(diamonds_by_cut, file = "optimization/data/diamonds_by_cut.csv")  
  
})




glimpse(mrgdf)

profvis({
a <- mrgdf %>%
  group_by(ftype) %>%
  summarise(nrec=n(), numret=sum(wt), nrecdiff=sum(taxbc!=taxbc.reform), nretdiff=sum(wt * (taxbc!=taxbc.reform)),
            taxbc=sum(taxbc * wt), taxbc.reform=sum(taxbc.reform * wt))
b <- a %>%
  mutate(diff=taxbc.reform - taxbc)
c <- b %>%
  write_csv("d:/tcdir/dan5.csv")
})


mrgdf %>%
  group_by(ftype) %>%
  summarise(nrec=n(), numret=sum(wt), nrecdiff=sum(taxbc!=taxbc.reform), nretdiff=sum(wt * (taxbc!=taxbc.reform)),
            taxbc=sum(taxbc * wt), taxbc.reform=sum(taxbc.reform * wt)) %>%
  mutate(diff=taxbc.reform - taxbc) %>%
  write_csv("d:/tcdir/dan1.csv")


mrgdf %>%
  select(ftype, RECID, wt, taxbc, taxbc.reform) %>%
  mutate(diff=taxbc.reform - taxbc) %>%
  write_csv("d:/tcdir/dan_recs.csv")

mrgdf %>%
  filter(taxbc != taxbc.reform) %>%
  mutate(diff=taxbc.reform - taxbc) %>%
  write_csv("d:/tcdir/dan_changerecs.csv")

mrgdf %>%
  group_by(ftype) %>%
  summarise(dupid=anyDuplicated(RECID))



# I had to add C:\Users\donbo\Anaconda3\Library\bin to the system path for Tax-Calculator to work!

# Don't change the environment variable
# path1 <- "C:/Users/donbo/Anaconda3"
# path2 <- "C:/Users/donbo/Anaconda3/Scripts"
# Sys.setenv(
#   PATH = paste(path1, path2,
#     Sys.getenv("PATH"), 
#     sep = ";"
#   )
# )
# 
# shell("echo %PATH% ", intern= TRUE)
tc.infile.fullpath <- shQuote(paste0(paste0(tc.dir, tc.fn)))

cmd1 <- "C:/Users/donbo/Anaconda3/Scripts/tc"
args <- c("D:/tcdir/synth10syn20.csv", "2013", 
          "--reform", "D:/Dropbox/RPrograms PC/OSPC/syndata4/tax_plans/brk4_1k_2013.json", 
          "--dump", 
          "--outdir", "D:/tcdir/")



infile <- "D:/tcdir/synth10syn20.csv"
taxplan <- shQuote(paste0(globals$taxplans.dir, reform, ".json"))
taxyear <- 2013

args <- c(infile, taxyear, 
          "--reform", taxplan,
          "--outdir", "D:/tcdir/")
system2(cmd1, args)



reform <- "brk4_1k_2013"
source.fn <- "synth10syn20.csv"

taxplan.fn <- paste0(reform, ".json")
# system.time(ref <- run_taxplan_getresults(source.fn, taxplan.fn))

# tc.wincmd <- function(tc.fn, tc.dir, tc.cli, taxyear=2013, taxplan.fn=NULL, taxplans.dir=NULL){
tc.dir <- globals$tc.dir
tc.fn <- source.fn
taxplans.dir <- globals$taxplans.dir

tc.cli <- globals$tc.cli
taxyear <- 2013

tc.infile.fullpath <- shQuote(paste0(paste0(tc.dir, tc.fn)))
tc.outdir <- shQuote(str_sub(tc.dir, 1, -2)) # must remove trailing "/" # djb CHANGE

taxplanstring <- NULL
if(!is.null(taxplan.fn)) taxplanstring <- paste0("--reform", " ", shQuote(paste0(paste0(taxplans.dir, taxplan.fn))))

cmd <- paste0(tc.cli, " ", tc.infile.fullpath, " ", taxyear, " ", taxplanstring, " ", "--dump --outdir ", tc.outdir)

# cmd <- tc.wincmd(source.fn, globals$tc.dir, globals$tc.cli, taxyear=2013, taxplan.fn=taxplan.fn, taxplans.dir=globals$taxplans.dir)

cmd <- shQuote("C:/Users/donbo/Anaconda3/Scripts/tc "D:/tcdir/synth10syn20.csv" 2013 --reform "D:/Dropbox/RPrograms PC/OSPC/syndata4/tax_plans/brk4_1k_2013.json" --dump --outdir "D:/tcdir/"")

# C:/Users/donbo/Anaconda3/Scripts/tc "D:/tcdir/synth10syn20.csv" 2013 --reform "D:/Dropbox/RPrograms PC/OSPC/syndata4/tax_plans/brk4_1k_2013.json" --dump --outdir "D:/tcdir"

p1 <- "C:/Users/donbo/Anaconda3/Scripts/tc"
p2 <- "D:/tcdir/synth10syn20.csv"
p3 <- 2013

p <- paste0(p1, " ", p2, " ", p3)
p

cmd2 <- paste0(tc.cli, " ", tc.infile.fullpath, " ", taxyear)
system(cmd2)


system(p)
system(shQuote(p))

ps <- shQuote(p, type="cmd2")
ps <- shQuote(p)
ps
system(ps)
shell(ps)

system("python --version")
system("where python")


cmd
system(cmd)

cmd1 <- "C:/Users/donbo/Anaconda3/Scripts/tc"
args <- c("D:/tcdir/synth10syn20.csv", "2013", 
          "--reform", "D:/Dropbox/RPrograms PC/OSPC/syndata4/tax_plans/brk4_1k_2013.json", 
          "--dump", 
          "--outdir", "D:/tcdir/")
system2(cmd1, args[1:2])


perlcmd <- 'print "Hello World\\n";'
## Not run: 
shell(shQuote(paste("perl -e", 
                    shQuote(perlcmd, type = "cmd")),
              type = "cmd2"))


system("dir")
system(shQuote("dir"))
shell("dir")



glimpse(mrgdf)

mrgdf %>%
  filter(c00100 > 250e3, taxbc>0) %>%
  group_by(ftype) %>%
  summarise(n=n(), 
            wt_sum=sum(wt), 
            affected_income_b=sum(10 * wt) / 1e9, 
            taxbc_b=sum(taxbc * wt) / 1e9,
            taxbc_reform_b=sum(taxbc.reform * wt) / 1e9) %>%
  mutate(tax_impact_b=taxbc_reform_b - taxbc_b,
         rough_overestimate_b=affected_income_b * .05) %>%
  as.data.frame


mrgdf %>%
  filter(c00100 > 250e3, taxbc>0) %>%
  group_by(ftype) %>%
  summarise(n=n(), 
            wt_sum=sum(wt),
            taxbc_b=sum(taxbc * wt) / 1e9, 
            affected_income_b=sum(100 * wt) / 1e9) %>%
  mutate(rough_overestimate_b=affected_income_b * .05) %>%
  as.data.frame


