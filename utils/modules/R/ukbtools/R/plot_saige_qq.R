
plot_saige_qq <- function(dt, 
                          aes_main = aes(x=-log10(pvalue.expected), y = -log10(pvalue.observed), color = csqs_category, label = NA),
                          aes_ribbon = aes(ymin=clower, ymax=cupper)){
    ggplot(dt, aes_main) +
            geom_ribbon(aes_ribbon, fill="grey80", color="grey80") + 
            geom_point_rast(size = 3) +
            geom_abline(linetype = 'dashed') + 
            geom_text_repel(color = 'black', size = 5) + 
            labs(title=phenotype, color=paste0("key_label")) +
            scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
            xlab(expression(paste(-log[10],'(Expected P-value)' ))) +
            ylab(expression(paste(-log[10],'(Observed P-value)' ))) +
            facet_wrap(~csqs_category) +
            theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
                  axis.title.y = element_text(margin=ggplot2::margin(r=10)),
                  plot.title = element_text(hjust=0.5))
    
}



