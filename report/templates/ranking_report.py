<%
        import os.path
        import warnings
        import pandas as pd
        import re
        import py_exp_calc.exp_calc as pxc
        import sys 
        sys.path.append("./report")
        import pyreport_helper as ph
        warnings.simplefilter(action='ignore', category=FutureWarning)

        number_of_positives = len(plotter.hash_vars["control_pos"])-1

        parse_name = {"string_ppi_combined": "STRING combined", 
                "string_ppi_textmining":"STRING textmining",
                "string_ppi_coexpression": "STRING coexpression",
                "string_ppi_neighborhood": "STRING neighborhood",
                "string_ppi_experimental": "STRING experiments",
                "string_ppi_cooccurence": "STRING cooccurrence",
                "phenotype": "HPO",
                "disease":"Disease",
                "pathway": "Pathway",
                "DepMap_effect_pearson": "DepMap Pearson",
                "string_ppi_database":"STRING databases",
                "DepMap_effect_spearman":"DepMap Spearman",
                "hippie_ppi": "Hippie",
                "DepMap_Kim":"DepMap Kim",
                "string_ppi_fusion":"STRING fusion",
                "gene_hgncGroup": "HGNC group",
                "integration_mean_by_presence": "IMP",
                "mean": "Mean",
                "max": "Max",
                "median": "Median",
                "el": "EL",
                "raw_sim": "GSM",
                "node2vec": "node2vec",
                "rf": "RF",
                "auc_down_ci_0.95": "AUROC-0.025",
                "auc_up_ci_0.95": "AUROC-0.975",
                "auc": "AUROC",
                } 

        def get_medianrank_size(var_name, groupby = ['annot_Embedding','annot','Embedding'], value = 'absolute_ranking', tops = [5,10,20]):
                df = pd.DataFrame(plotter.hash_vars[var_name][1:], columns = plotter.hash_vars[var_name][0])
                df_final = df.groupby(groupby)[value].size().reset_index()
                for top in tops:
                        df_top = df.groupby(groupby)[value].apply(lambda x: sum(pd.to_numeric(x)<=top)).reset_index()
                        df_top.rename(columns={value: f"Top {top}"}, inplace=True)
                        df_final = pd.concat([df_final, df_top[[f"Top {top}"]]], axis=1)
                #col_names = plotter.hash_vars[var_name][0]
                #col_names.append("size")
                return [df_final.columns.tolist()] + df_final.values.tolist()

        for table in plotter.hash_vars.keys():
                if table == "parsed_non_integrated_rank_summary" or table == "parsed_integrated_rank_summary":
                        ph.parse_table(plotter, table, parse_name, include_header=True)
                else:
                        ph.parse_table(plotter, table, parse_name)
        if plotter.hash_vars.get('parsed_non_integrated_rank_summary') is not None:
                ph.order_columns(plotter,'parsed_non_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_integrated_rank_summary') is not None:
                ph.order_columns(plotter,'parsed_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_non_integrated_rank_pos_cov') is not None:
                ph.order_columns(plotter,'parsed_non_integrated_rank_pos_cov',0)
                plotter.hash_vars["parsed_non_integrated_rank_pos_cov"] = ph.modify_by_cols(plotter,"parsed_non_integrated_rank_pos_cov", [3], lambda x: float(x)/number_of_positives * 100)
                plotter.hash_vars['parsed_non_integrated_rank_pos_cov'] = ph.parse_heatmap_from_flat(plotter.hash_vars['parsed_non_integrated_rank_pos_cov'][1:],1,2,3,None,None)

        if plotter.hash_vars.get('parsed_integrated_rank_pos_cov') is not None:
                ph.order_columns(plotter,'parsed_integrated_rank_pos_cov',0)
                plotter.hash_vars["parsed_integrated_rank_pos_cov"] = ph.modify_by_cols(plotter,"parsed_integrated_rank_pos_cov", [3], lambda x: float(x)/number_of_positives * 100)
                plotter.hash_vars['parsed_integrated_rank_pos_cov'] = ph.parse_heatmap_from_flat(plotter.hash_vars['parsed_integrated_rank_pos_cov'][1:],1,2,3,None,None)

        if plotter.hash_vars.get('parsed_annotation_grade_metrics') is not None:
                ph.order_columns(plotter,'parsed_annotation_grade_metrics',0)

        if plotter.hash_vars.get("non_integrated_rank_cdf"):
                plotter.hash_vars["non_integrated_tops"] = get_medianrank_size("non_integrated_rank_cdf", tops=[5,10,20,100])
        
        if plotter.hash_vars.get("integrated_rank_cdf"):
                plotter.hash_vars["integrated_tops"] = get_medianrank_size("integrated_rank_cdf", ["integration_Embedding","integration","Embedding"], tops = [5,10,20,100])
%>

<% plotter.set_header() %>

<% txt="From eGSM to Backup Genes." %>
${plotter.create_title(txt, id='main_backup_gene', hlevel=1, indexable=True, clickable=False)}


<% txt="Workflow for Backup Gene Benchmarking" %>
${plotter.create_title(txt, id='main_backup_gene', hlevel=2, indexable=True, clickable=False)}

<%
        graph=f"""
        graph LR
            SI[Manually curated <br>from literature]
            Papi[<span style="color:#280054">Big Papi </span>]
            Digenic[<span style="color:#280054">Digenic </span>]
            P[<span style="color:#023020">Positive</span>]
            N[<span style="color:color:#500000">Negative</span>]
            CDF[<span style="color:#000000">CDF</span>]
            ROC[<span style="color:#000000">ROC </span>]
            Coverage[<span style="color:#000000">Coverage</span>]
            subgraph C [<u><b>Backup <br> GoldStandard</b></u>]
            P
            N
            end
            subgraph SII [Double knock-out screenings]
            Papi
            Digenic
            end
            SI--> P
            SII -- Significative interactions \\n p-value cutoff .05 --> C
            P--> CDF
            P--> Coverage
            P & N --> ROC
            style P fill:#C0EAB9,stroke:#0D3E05,stroke-width:3px
            style N fill:#FF8C80,stroke:#8C1F14,stroke-width:3px
            style CDF fill:#D4D4D4,stroke:#7B7B7B,stroke-width:2px
            style ROC fill:#D4D4D4,stroke:#7B7B7B,stroke-width:2px
            style Coverage fill:#D4D4D4,stroke:#7B7B7B,stroke-width:2px
            style Papi fill:#547DC6,stroke:#333,stroke-width:2px
            style Digenic fill:#547DC6,stroke:#333,stroke-width:2px
            style SII fill:#ABC9FF,stroke:#333,stroke-width:2px
            style SI fill:#ABC9FF,stroke:#333,stroke-width:2px
            style C fill:#FFFCDE,stroke:#333,stroke-width:2px,text-align:center
        """
%>
${plotter.mermaid_chart(graph)}
${ph.make_title(plotter,"figure", "seed_wflow", """Workflow of the benchmarking process. The reference for gene backup controls have been 
        obtained either from the scientific literature or from double CRIPSR knock-out screenings. The last ones were selected based on p-value on most tested cellular lines. 
        """)} 

<% txt="Backup Gene Coverage" %>
${plotter.create_title(txt, id='backup_cov', hlevel=2, indexable=True, clickable=False)}

<div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
        <div style="margin-right: 10px;">
                % if plotter.hash_vars.get('parsed_non_integrated_rank_pos_cov') is not None:
                        ${ plotter.heatmap(id = 'parsed_non_integrated_rank_pos_cov', title="",header = True, row_names = True, 
                                config= {"setMinX":0,
                                "setMaxX":100, 
                                "xAxisTitle": "Coverage", 
                                "samplesClustered":True,
                                "showSmpDendrogram":False}) }
                % endif
        </div>
        <div style="margin-left: 10px;"> 
                % if plotter.hash_vars.get('parsed_integrated_rank_pos_cov') is not None: 
                        ${ plotter.heatmap(id = 'parsed_integrated_rank_pos_cov', header = True, title="", row_names = True, config= {"setMinX":0,"setMaxX":100, "xAxisTitle": "Coverage", "samplesClustered":True,"showSmpDendrogram":False}) }
                % endif
        </div>
</div>
${ph.make_title(plotter,"figure", "coverage_bars", """Coverage obtained in each individual (A)
 or integrated (B) eGSM. In both plots, x axis reflects the number of positive control genes with information on the adjacency matrix, discarding does 
 with zero or minimum value on edges for the corresponding seed.""")}

<% txt="Performance metrics" %>
${plotter.create_title(txt, id='perf_metrics', hlevel=2, indexable=True, clickable=False)}

<% txt="Rank distributions" %>
${plotter.create_title(txt, id='summ_rank_dis', hlevel=3, indexable=True, clickable=False)}

<div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
        <div style="margin-right: 10px;">
                        % if plotter.hash_vars.get('non_integrated_rank_cdf') is not None: 
                                ${plotter.boxplot(id= 'non_integrated_rank_cdf', header= True, row_names= False, default= False, fields= [5],  smp_attr= [0,1,2], group = "Embedding",
                                   title= "(A) Individual eGSM",
                                        x_label= "Normalized rank",
                                        config= {
                                                "graphOrientation": "vertical",
                                                "colorBy" : "Embedding",
                                                "groupingFactors" :
                                                ["Embedding"],
                                                "titleFontStyle": "italic",
                                                "titleScaleFontFactor": 0.7,
                                                "smpTextRotate": 45,
                                                "segregateSamplesBy": "annot"})}
                        % endif
        </div>
        <div style="margin-left: 10px;">
                        % if plotter.hash_vars.get('integrated_rank_cdf') is not None: 
                                ${plotter.boxplot(id= 'integrated_rank_cdf', header= True, row_names= False, default= False, fields = [5], smp_attr= [0,1,2], group= "Embedding", 
                                        title= "(B) Integrated eGSM",
                                        xlabel= "Normalized rank",
                                        config= {
                                                "graphOrientation": "vertical",
                                                "colorBy" : "Embedding",
                                                "xAxisTitle": "Normalized rank",
                                                "groupingFactors" :
                                                ["Embedding"],
                                                "titleFontStyle": "italic",
                                                "titleScaleFontFactor": 0.7,
                                                "smpTextRotate": 45,
                                                "segregateSamplesBy": "integration"})}
                        % endif
        </div>
</div>
${ph.make_title(plotter,"figure", "rank_boxplot", f"""Rank distributions in each individual (A)
 or integrated (B) eGSM. In both plots, y axis ({ph.italic("Normalized ranks")}) represent the rank is normalized on 0-1 range.""")}

<% txt="Tops" %>
${plotter.create_title(txt, id='tops', hlevel=3, indexable=True, clickable=False)}
<%text=[]%>
<% text.append(plotter.barplot(id="non_integrated_tops", header=True, fields=[2,4,5,6,7], smp_attr=[1,2], 
        title= "(A) Individual eGSM", 
        x_label="Number of top-ranked Backups",
        config={
        'graphOrientation' : 'vertical',
        'segregateSamplesBy': "annot",
        "titleFontStyle": "italic",
        'setMaxX': 87,
        "smpTextRotate": 45,
        "fontScaleFontFactor": 1.5,
        "axisTitleScaleFontFactor": 1.5,
        "axisTickScaleFontFactor": 1.5,
        "legendTextScaleFontFactor": 1.2,
        "legendTitleScaleFontFactor": 1.4
        })) %>
<% text.append(plotter.barplot(id="integrated_tops", header=True, fields=[2,4,5,6,7], smp_attr=[1,2], 
        title= "(B) Integrated eGSM",
        x_label="Number of top-ranked Backups", 
        config={
        'graphOrientation' : 'vertical',
        'segregateSamplesBy': "integration",
        "titleFontStyle": "italic",
        "titleScaleFontFactor": 0.7,
        "smpTextRotate": 45,
        "smpLabelScaleFontFactor": 0.3,
        'setMaxX': 87,
        "fontScaleFontFactor": 1.5,
        "axisTitleScaleFontFactor": 1.5,
        "axisTickScaleFontFactor": 1.5,
        "legendTextScaleFontFactor": 1.2,
        "legendTitleScaleFontFactor": 1.4
        })) %>
<% text.append(ph.make_title(plotter,"figure", "agg_tops", f"""Top 5,10,20,100 on different individual (A) and integrated (B) eGSM."""))%>
${ph.collapsable_data(plotter,"Tops absolute values", None, "tops_absolute", "\n".join(text))}

<div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
        <%  

        plotter.hash_vars["non_integrated_tops_relative"] = plotter.hash_vars["non_integrated_tops"][:]
        for idx, row in enumerate(plotter.hash_vars["non_integrated_tops"][1:]):
                for col in [4,5,6,7]: 
                        plotter.hash_vars["non_integrated_tops_relative"][idx+1][col] = row[col]/number_of_positives * 100

        plotter.hash_vars["integrated_tops_relative"] = plotter.hash_vars["integrated_tops"][:]
        for idx, row in enumerate(plotter.hash_vars["integrated_tops"][1:]):
                for col in [4,5,6,7]:
                        plotter.hash_vars["integrated_tops_relative"][idx+1][col] = row[col]/number_of_positives * 100
        %>
        ${plotter.barplot(id="non_integrated_tops_relative", header=True, fields=[2,4,5,6,7], smp_attr=[1,2], 
                title= "(A) Individual eGSM", 
                x_label=f"% of top-ranked Backups",
                config={
                'graphOrientation' : 'vertical',
                'segregateSamplesBy': "annot",
                "titleFontStyle": "italic",
                'setMaxX': 100,
                "smpTextRotate": 45,
                "fontScaleFontFactor": 1.5,
                "axisTitleScaleFontFactor": 1.5,
                "axisTickScaleFontFactor": 1.5,
                "legendTextScaleFontFactor": 1.2,
                "legendTitleScaleFontFactor": 1.4
                })}
        ${plotter.barplot(id="integrated_tops_relative", header=True, fields=[2,4,5,6,7], smp_attr=[1,2], 
                title= "(B) Integrated eGSM",
                x_label=f"% of top-ranked Backups", 
                config={
                'graphOrientation' : 'vertical',
                'segregateSamplesBy': "integration",
                "titleFontStyle": "italic",
                "titleScaleFontFactor": 0.7,
                "smpTextRotate": 45,
                'setMaxX': 100,
                "fontScaleFontFactor": 1.5,
                "axisTitleScaleFontFactor": 1.5,
                "axisTickScaleFontFactor": 1.5,
                "legendTextScaleFontFactor": 1.2,
                "legendTitleScaleFontFactor": 1.4
                })}
</div>
${ph.make_title(plotter,"figure", "agg_tops", f"""Top 5,10,20,100 on different individual (A) and integrated (B) eGSM.""")}
<% txt="Performance curves" %>
${plotter.create_title(txt, id='perf_curves', hlevel=2, indexable=True, clickable=False)}

<% txt="CDF" %>
${plotter.create_title(txt, id='cdf_curves', hlevel=3, indexable=True, clickable=False)}

<a href="https://academic.oup.com/bib/article/23/2/bbac019/6521702#330302198">Xiao Yuan et al. Evaluation of phenotype-driven gene prioritization methods for Mendelian diseases, Briefings in Bioinformatics, Volume 23, Issue 2, March 2022, bbac019 </a>
<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="non_integrated_rank_cdf", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                width=600, height=600, matplot_width = 6, matplot_height = 6, dpi=1200,  plotting_function= lambda data, plotter_list: ph.plot_with_facet(plot_type="ecdf",data=data, 
                                        plotter_list=plotter_list, x="rank", col="annot", 
                                        hue="Embedding", col_wrap=3, 
                                        suptitle="", x_label="Normalized Rank", y_label="TPR", top=0.9))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="integrated_rank_cdf", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                width=600, height=600, matplot_width = 6, matplot_height = 6, dpi=1200, plotting_function= lambda data, plotter_list: ph.plot_with_facet(plot_type="ecdf",data=data, 
                                        plotter_list=plotter_list, x="rank", 
                                        col="integration", hue="Embedding", col_wrap=2, suptitle="", x_label="Normalized Rank", y_label="TPR", top=0.8))}
        % endif
</div>
${ph.make_title(plotter,"figure", "cdf_curve", f"""CDF curves by each individual (A)
 or integrated (B) eGSM. In both plots, y axis represent 
 the true positive rate ({ph.italic("TPR")}) and x axis ({ph.italic("Normalized Rank")}) the rank normalized from 0 to 1.""")}

<% txt="ROC" %>
${plotter.create_title(txt, id='roc_curves', hlevel=3, indexable=True, clickable=False)}

<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_measures") is not None: 

                 ${ plotter.static_plot_main( id="non_integrated_rank_measures", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                width=600, height=600, matplot_width = 6, matplot_height = 6, dpi=1200,plotting_function= lambda data, plotter_list: ph.plot_with_facet(plot_type="lineplot", data=data,
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='annot', 
                                        hue='Embedding', col_wrap=3, suptitle="", 
                                        top=0.9, x_label="FPR", y_label="TPR"))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_measures") is not None: 
                 ${ plotter.static_plot_main( id="integrated_rank_measures", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6], 
                                width=600, height=600, matplot_width = 6, matplot_height = 6, dpi=1200,plotting_function= lambda data, plotter_list: ph.plot_with_facet(plot_type="lineplot",data=data, 
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='integration', 
                                        hue='Embedding', col_wrap=2, suptitle="", 
                                        top=0.8, labels = 'Embedding', x_label="FPR", y_label="TPR"))}
        % endif
</div>
${ph.make_title(plotter,"figure", "roc_curve", f"""ROC in each individual (A) or integrated (B) eGSM.""")}

<%txt=[]%>
% if plotter.hash_vars.get("parsed_non_integrated_rank_summary") is not None: 
        <% txt.append(plotter.line(id= "parsed_non_integrated_rank_summary", fields= [1, 7, 13, 8], header= True, row_names= True, smp_attr=[0,2],
                responsive= False,
                height= '400px', width= '400px', x_label= 'AUROC',
                title= "(A) Individual eGSM",
                config= {
                        'showLegend' : True,
                        'graphOrientation' : 'vertical',
                        "titleFontStyle": "italic",
                        "titleScaleFontFactor": 0.7,
                        'setMinX': 0,
                        'setMaxX': 1,
                        "smpTextRotate": 45,
                        "segregateSamplesBy": "Embedding"
                        })) %>
% endif
% if plotter.hash_vars.get('parsed_integrated_rank_summary') is not None: 
        <% txt.append(plotter.line(id= "parsed_integrated_rank_summary", fields=  [2, 7, 13, 8], header= True, row_names= True, smp_attr = [0,1],
                responsive= False,
                height= '400px', width= '400px', x_label= 'AUROC',
                title= "(B) Integrated eGSM",
                config= {
                        'showLegend' : True,
                        'graphOrientation' : 'vertical',
                        "titleFontStyle": "italic",
                        "titleScaleFontFactor": 0.7,
                        'setMinX': 0.5,
                        'setMaxX': 1,
                        'stripTextScaleFontFactor': 1.7,
                        'xAxisTextScaleFontFactor': 1.7,
                        'xAxisTitleScaleFontFactor': 1.8,
                        'legendTextScaleFontFactor': 1.4,
                        "smpTextRotate": 45,
                        "segregateSamplesBy": "Integration"
                        })) %>
% endif
<% txt.append(ph.make_title(plotter,"figure", "roc_ic", f"""AUROC Confidence Interval (CI) in each individual (A) or integrated (B) eGSM. CI was obtained by a 1000 iteration bootstrap.""")) %>
${ph.collapsable_data(plotter,"AUROC Confidence Interval", None, "auroc_ci", "\n".join(txt))}