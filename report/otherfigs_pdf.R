#' ---
#' title: Data preparation, Global TB Report 2022
#' author: Philippe Glaziou (updated by Irwin in 2022)
#' date: 10/08/2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' # Preamble
#' (Last updated: `r Sys.Date()`)
#'
#' Export files to GTB database
#'
library(data.table)
library(here)
library(gtbreport)
library(whomap)

source(here('report/ch1_2_dataprep.R'))


#GTB report 2022
#FIG.18 and FIG.19

#COUNTRIES WHICH MET THE 2020 MILESTONE IN 2021 FOR TB DEATHS

p18_data<- est %>%
  filter (country== "Bangladesh"
          |country== "Ethiopia" 
          |country== "Kenya"
          |country== "Mozambique"
          |country== "Uganda"
          |country== "United Republic of Tanzania"
          |country== "Zambia"
          |country== "Russian Federation")

p18<-qplot(year,
           mort.num/ 1e3,
           data = p18_data,
           geom = 'line',
           colour = I('black')) +
  geom_ribbon(aes(year, ymin = mort.lo.num / 1e3, ymax = mort.hi.num / 1e3),
              fill = I('orange'),
              #fill=gtbreport::palette_gtb("inc"),
              alpha = 0.4) +
  geom_hline(aes(yintercept = mort.milestone/1e3), linetype = I(2)) +
  facet_wrap(~country, nrow=2, scales="free",
             # Use the labeller function to make sure long country names are wrapped in panel headers
             labeller = label_wrap_gen(width = 25)) +
  #  expand_limits(y = 1) +
  xlab('') + ylab('TB deaths (total, in thousands) per year (log scale)') +
  scale_y_log10() +
  theme_gtb()

p18

write.csv(p18_data[, .(year, iso3, country, mort.num, mort.lo.num, mort.hi.num, mort.milestone)], 
          file=here('report/pdf/fig18.csv'))

ggsave(here('report/pdf/fig18.pdf'), width=12, height=8)


#COUNTRIES WHICH MET THE 2020 MILESTONE IN 2021 FOR TB INCIDENCE RATE

p19_data<- est %>%
  filter(country== "Ethiopia" 
         |country== "Kenya"
         |country== "Lesotho"
         |country== "Namibia"
         |country== "South Africa"
         |country== "United Republic of Tanzania"
         |country== "Zambia"|country== "Cambodia"
         |country== "Russian Federation"
         |country== "Zimbabwe")

p19<-qplot(year,
           inc,
           data = p19_data,
           geom = 'line',
           colour = I('black')) +
  geom_ribbon(aes(year, ymin = inc.lo, ymax = inc.hi),
              #fill = I('seagreen'),
              fill=gtbreport::palette_gtb("inc"),
              alpha = 0.4) +
  #geom_line(aes(year, newinc)) +
  geom_hline(aes(yintercept = inc.milestone), linetype = I(2)) +
  #facet_wrap(~ country, scales = 'free_y', ncol = 5) +
  facet_wrap(~country, nrow=3, scales="free",
             # Use the labeller function to make sure long country names are wrapped in panel headers
             labeller = label_wrap_gen(width = 25)) +
  #  expand_limits(y = 1) +
  xlab('') + ylab('Rate per 100 000 population per year (log scale)') +
  scale_y_log10() +
  theme_gtb()

p19

write.csv(p19_data[, .(year, iso3, country, inc, inc.lo, inc.hi, inc.milestone)], 
          file=here('report/pdf/fig19.csv'))

ggsave(here('report/pdf/fig19.pdf'), width=12, height=8)





#2021 REPORT CODE BY PHILIPPE BELOW


# fig 17
hest[year>=2015, dec.inc := 1-inc/inc[1], by=iso3]
lst <- hest[year==2021, .(iso3, dec.inc)][order(-dec.inc)][1:10, iso3]
sel <- hest$iso3 %in% lst & hest$year>2014

qplot(
  year,
  inc,
  data = hest[sel],
  geom = 'line',
  colour = I('grey20')
) +
  geom_ribbon(
    aes(year, ymin = inc.lo, ymax = inc.hi),
    fill = I('blue'),
    alpha = 0.4
  ) +
  geom_hline(aes(yintercept = inc.milestone), linetype = I(2)) +
  facet_wrap(~country, nrow=2, scales="free",
             labeller = label_wrap_gen(width = 25)) +
  xlab('') + ylab('TB incidence rate per 100 000 per year (log scale)') +
  scale_y_log10() +
  theme_gtb()

ggsave(here('report/pdf/fig17.pdf'), width=12, height=8)
fwrite(hest[sel, .(iso3,year,inc,inc.lo,inc.hi,inc.milestone)],
       file=here('report/pdf/fig17.csv'))




# fig 18
hest[year>=2015, dec.mort := 1-mort.num/mort.num[1], by=iso3]
lst <- hest[year==2021, .(iso3, dec.mort)][order(-dec.mort)][1:10, iso3]
sel <- hest$iso3 %in% lst & hest$year>2014

qplot(
  year,
  mort.num / 1e3,
  data = hest[sel],
  geom = 'line',
  colour = I('grey20')
) +
  geom_ribbon(
    aes(year, ymin = mort.lo.num / 1e3, ymax = mort.hi.num / 1e3),
    fill = I('blue'),
    alpha = 0.4
  ) +
  geom_hline(aes(yintercept = mort.milestone / 1e3), linetype = I(2)) +
  facet_wrap(~country, nrow=2, scales="free",
             labeller = label_wrap_gen(width = 25)) +
  xlab('') + ylab('TB deaths (total, in thousands) per year (log scale)') +
  scale_y_log10() +
  theme_gtb()

ggsave(here('report/pdf/fig18.pdf'), width=12, height=8)
fwrite(hest[sel, .(iso3,year,mort.num,mort.lo.num,mort.hi.num,mort.milestone)],
       file=here('report/pdf/fig18.csv'))




# fig 20
qplot(
    year,
    inc.num / 1e6,
    data = global,
    geom = 'line',
    colour = I('blue')
  ) +
  geom_ribbon(
    aes(year, ymin = inc.lo.num / 1e6, ymax = inc.hi.num / 1e6),
    fill = I('blue'),
    alpha = 0.4
  ) +
  geom_line(aes(year, c.newinc / 1e6)) +
  ylab('Millions per year (log scale)') + xlab('') +
  expand_limits(y = 1) +
  scale_y_log10() +
  theme_bw(base_size = 14)

ggsave(here('report/pdf/fig20.pdf'), width=8.3, height=5.8)
fwrite(global[, .(year,inc.num,inc.lo.num,inc.hi.num,c.newinc)],
       file=here('report/pdf/fig20.csv'))








