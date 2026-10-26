# SMPM E4.5 S1 — corrected publication figures (no re-estimation)
# Run from the repository root after running the main application.
required <- c('ggplot2','dplyr','tidyr','readr','gridExtra')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install packages: install.packages(c(',paste(sprintf('"%s"',missing),collapse=', '),'))')
library(ggplot2); library(dplyr); library(tidyr); library(readr)
root <- getwd()
find_input <- function(name, folders) {
  paths <- file.path(root, folders, name)
  paths <- c(paths, file.path(root,name))
  hits <- paths[file.exists(paths)]
  if(!length(hits)) stop('Missing file: ',name,'; searched: ',paste(paths,collapse=' | '))
  hits[1]
}
score_file <- find_input('cell_scores.csv',c('results'))
eigen_file <- find_input('supervised_eigenfunctions.csv',c('results'))
out <- file.path(root,'figures'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
validation_out <- file.path(root,'validation'); dir.create(validation_out,recursive=TRUE,showWarnings=FALSE)

scores <- read_csv(score_file,show_col_types=FALSE)
stopifnot(all(c('celltype','SMPC1','SMPC2') %in% names(scores)), nrow(scores)==2015L)
stopifnot(!anyNA(scores[,c('celltype','SMPC1','SMPC2')]))
scores$celltype <- factor(scores$celltype,levels=c('EPI','PE','TE'))
stopifnot(!anyNA(scores$celltype))
centers <- scores %>% group_by(celltype) %>% summarise(SMPC1=mean(SMPC1),SMPC2=mean(SMPC2),n=n(),.groups='drop')
write_csv(centers,file.path(validation_out,'Figure1_centroids.csv'))
colors <- c(EPI='#2166AC',PE='#D17B27',TE='#238B75')
# Filled shapes allow a black outline for the mean markers.
shapes <- c(EPI=21,PE=24,TE=22)
zoom_x <- as.numeric(quantile(scores$SMPC1,c(.02,.98),na.rm=TRUE))
zoom_y <- as.numeric(quantile(scores$SMPC2,c(.02,.98),na.rm=TRUE))
# Always include all three means inside the zoom window.
zoom_x <- range(c(zoom_x,centers$SMPC1)); zoom_y <- range(c(zoom_y,centers$SMPC2))
pad <- function(v) { d <- diff(v); v+c(-.08,.08)*d }
zoom_x <- pad(zoom_x); zoom_y <- pad(zoom_y)

# Individual cells: solid points with mapped colour (NOT fill).
# Shape 16 does not use the fill aesthetic; using fill with shape 16
# was why the 2,015 observations disappeared in the previous PDF.
base_scores <- ggplot(scores, aes(SMPC1, SMPC2)) +
  geom_point(aes(color=celltype), shape=16, alpha=.22, size=1.15,
             show.legend=TRUE) +
  stat_ellipse(aes(group=celltype,color=celltype), type='norm', level=.68,
               linewidth=.7, show.legend=FALSE) +
  geom_point(data=centers,
             aes(SMPC1, SMPC2, fill=celltype, shape=celltype),
             inherit.aes=FALSE, color='black', stroke=1.0, size=4.2,
             show.legend=FALSE) +
  scale_color_manual(values=colors, name='Lineage', drop=FALSE) +
  scale_fill_manual(values=colors, guide='none', drop=FALSE) +
  scale_shape_manual(values=shapes, guide='none', drop=FALSE) +
  labs(x='SMPM 1 score',y='SMPM 2 score') +
  theme_classic(base_size=11) +
  theme(legend.position='top',legend.title=element_blank(),
        plot.title=element_text(face='bold',size=11),
        axis.title=element_text(face='bold')) +
  guides(color=guide_legend(override.aes=list(alpha=1,size=3)))
full <- base_scores + ggtitle('(A) All cells: full range')
zoom <- base_scores + coord_cartesian(xlim=zoom_x,ylim=zoom_y,expand=FALSE) +
  ggtitle('(B) Central region')
# coord_cartesian zoom retains the entire dataset for ellipses and means.
combined <- gridExtra::arrangeGrob(full,zoom,ncol=2)

# Figure 2: compute true grid differences, then aggregate by the right endpoint
# into 1-Mb display windows. This is a display summary, not peak calling.
eigen <- read_csv(eigen_file,show_col_types=FALSE)
marks <- c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3')
stopifnot(all(c('mark','position_bp','SMPC1','SMPC2') %in% names(eigen)))
stopifnot(setequal(unique(eigen$mark),marks),all(table(eigen$mark)==1001L))
long <- eigen %>% pivot_longer(c(SMPC1,SMPC2),names_to='component',values_to='phi') %>%
  mutate(component=factor(component,levels=c('SMPC1','SMPC2'),labels=c('SMPM 1','SMPM 2')),
         mark=factor(mark,levels=marks)) %>%
  arrange(component,mark,position_bp) %>% group_by(component,mark) %>%
  mutate(left_bp=lag(position_bp),signed_increment=phi-lag(phi)) %>%
  ungroup() %>% filter(!is.na(left_bp))
stopifnot(nrow(long)==10000L,all(long$position_bp>long$left_bp))
write_csv(long %>% select(component,mark,left_bp,position_bp,signed_increment),
          file.path(validation_out,'Figure2_exact_grid_increments.csv'))
region <- long %>% mutate(window=pmax(0,ceiling(position_bp/1e6)-1L),
                          midpoint_mb=window+0.5) %>%
  group_by(component,mark,window,midpoint_mb) %>%
  summarise(signed_mass=sum(signed_increment),.groups='drop')
write_csv(region,file.path(validation_out,'Figure2_signed_increments_1Mb.csv'))
check <- long %>% group_by(component,mark) %>% summarise(grid_total=sum(signed_increment),.groups='drop') %>%
  left_join(region %>% group_by(component,mark) %>% summarise(window_total=sum(signed_mass),.groups='drop'),
            by=c('component','mark'))
stopifnot(max(abs(check$grid_total-check$window_total))<1e-9)
write_csv(check,file.path(validation_out,'Figure2_mass_check.csv'))

# facet_grid(scales='free_y') uses a common vertical scale ACROSS BOTH
# components of the same histone mark (one scale per row).
fig2 <- ggplot(region,aes(midpoint_mb,signed_mass)) +
  geom_hline(yintercept=0,linewidth=.35,color='#404040') +
  geom_col(aes(fill=signed_mass>=0),width=.94,show.legend=FALSE) +
  scale_fill_manual(values=c('TRUE'='#2478A3','FALSE'='#BB5D4B')) +
  facet_grid(mark~component,scales='free_y') +
  scale_x_continuous(breaks=seq(0,200,25),expand=expansion(mult=c(0,.01))) +
  labs(x='Chromosome 1 position (Mb)',
       y='Signed principal-measure increment (1-Mb windows)') +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill='#F0F0F0',color='#BDBDBD'),
        strip.text=element_text(face='bold'),axis.title=element_text(face='bold'))

save_both <- function(plot,stem,w,h) {
  ggsave(file.path(out,paste0(stem,'.pdf')),plot,width=w,height=h,units='in',device='pdf',bg='white')
  if(requireNamespace('ragg',quietly=TRUE)) {
    ggsave(file.path(out,paste0(stem,'.png')),plot,width=w,height=h,units='in',dpi=400,device=ragg::agg_png,bg='white')
  } else {
    ggsave(file.path(out,paste0(stem,'.png')),plot,width=w,height=h,units='in',dpi=400,bg='white')
  }
}
save_both(combined,'Figure2_SMPM_scores',11.5,5.0)
save_both(fig2,'Figure3_SMPM_signed_increments',9.0,9.2)
cat('Figures written to: ',normalizePath(out),'\n',sep='')
