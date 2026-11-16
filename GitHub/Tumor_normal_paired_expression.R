library(ggplot2)
library(ggpubr)
#1
df <- read.table("box2.txt", header = T)
ggplot(data=df, aes(x = Group, y = Expression,
                    color=Group)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  scale_color_manual(limits=c("Normal","Tumor"),
                     values=c("#E29827","#922927"))+
  stat_compare_means(method = "t.test",paired = TRUE,
                     comparisons=list(c("Normal", "Tumor")))+
  geom_jitter(alpha = 0.3,size=3, aes(fill=Patient),shape=21)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  facet_wrap(~Gene, scales = "free_y")+
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')+
  labs(y='Expression')

ggplot(data=df, aes(x = Group, y = Expression,
                    color=Group)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  scale_color_manual(limits=c("Normal","Tumor"),
                     values=c("#4E79A7","#E15759"))+
  stat_compare_means(method = "t.test",paired = TRUE,
                     comparisons=list(c("Normal", "Tumor")))+
  geom_jitter(alpha = 0.3,size=3, aes(fill=Patient),shape=21)+
  geom_line(aes(group = Patient),
            color = 'grey40', lwd = 0.5)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  facet_wrap(~Gene, scales = "free_y")+
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')+
  labs(y='Expression')

p<-ggplot(data=df, aes(x = Group, y = Expression,
                    color=Group)) +
  geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA)+
  scale_color_manual(limits=c("Normal","Tumor"),
                     values=c("#4E79A7","#E15759"))+
  stat_compare_means(method = "t.test",paired = TRUE,
                     comparisons=list(c("Normal", "Tumor")))+
  geom_jitter(size=1.5, aes(fill=Patient),
              shape=21,position = position_dodge(0.5))+
  geom_line(aes(group = Patient),
            color = '#A5AA99', lwd = 0.5,position = position_dodge(0.5))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  facet_wrap(~Gene, scales = "free_y")+
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'right')+
  labs(y='Expression')
ggsave(p, file="box2.pdf",width=4.5,height=235)
