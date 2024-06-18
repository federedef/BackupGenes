<%
        import os.path
        import warnings
        import pandas as pd
        import re
        import py_exp_calc.exp_calc as pxc
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
                "auc_down_ci_0.95": "AUROC-down",
                "auc_up_ci_0.95": "AUROC-up",
                "auc": "AUROC",
                }

        # Text
        #######

        def parse_heatmap_from_flat(data,nrow,ncol,nvalue):
                pairs = {}
                for row in data:
                        if not pairs.get(row[nrow]):
                                pairs[row[nrow]] = {}
                        pairs[row[nrow]][row[ncol]] = row[nvalue]
                mat, row, col = pxc.to_wmatrix_rectangular(pairs)
                table = [["-",*col]]
                for idx,elem in enumerate(row):
                        table.append([elem, *mat[idx,:].tolist()])
                return table

        def italic(txt):
                return f"<i>{txt}</i>"

        def collapsable_data(click_title, click_id, container_id, txt, indexable=False, hlevel=1):
                collapsable_txt = f"""
                {plotter.create_title(click_title, id=click_id, indexable=indexable, clickable=True, hlevel=hlevel, t_id=container_id)}\n
                <div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
                        {plotter.create_collapsable_container(container_id, txt)}
                </div>"""
                return collapsable_txt

        def make_title(type, id, sentence):
                if type == "table":
                        key = f"tab:{id}"
                        html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {plotter.add_table(key)} </b> {sentence} </p>"
                elif type == "figure":
                        key = id
                        html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {plotter.add_figure(key)} </b> {sentence} </p>"
                return html_title

        def parsed_string(data, blacklist = ["sim"]):
                words = []
                for word in data.split("_"):
                        for blackword in blacklist:
                                word = re.sub(blackword,"",word)
                        word = word.capitalize()
                        words.append(word)
                parsed_data = " ".join(words)
                return parsed_data

        def parse_data(table, blacklist = ["sim"], column = "all"):
                parsed_table = []
                for i,row in enumerate(table):
                        parsed_table.append(row)
                        for j,data in enumerate(row):
                                if type(data) == str and not data.startswith("HGNC:"):
                                        if parse_name.get(data):
                                                parsed_table[i][j] = parse_name[data]
                                        else:
                                                parsed_table[i][j] = parsed_string(data, blacklist)
                                else:
                                        continue
                return parsed_table
                
        def order_columns(name, column):
                tab_header = plotter.hash_vars[name].pop(0)
                plotter.hash_vars[name].sort(key=lambda x: x[column])
                plotter.hash_vars[name].insert(0, tab_header)

        def parse_table(name, blacklist=["sim"], include_header = False):
                if not include_header:
                        tab_header = plotter.hash_vars[name].pop(0)
                        plotter.hash_vars[name] = parse_data(plotter.hash_vars[name])
                        if name == "parsed_integrated_rank_summary": print(plotter.hash_vars[name])
                        plotter.hash_vars[name].insert(0, tab_header)
                else:
                        plotter.hash_vars[name] = parse_data(plotter.hash_vars[name])

        def plot_with_facet(data, plotter_list, plot_type="", x='fpr', y='tpr', col=None, hue=None, col_wrap=4, suptitle=None, top=0.7, labels = None, x_label=None, y_label=None):
                if plot_type == "scatterplot":
                        g = plotter_list["sns"].FacetGrid(data, col_wrap=col_wrap, col=col, hue=hue, aspect=1).map(plotter_list["sns"].scatterplot, x, y)
                elif plot_type == "lineplot":
                        g = plotter_list["sns"].FacetGrid(data, col_wrap=col_wrap, col=col, hue=hue, aspect=1).map(plotter_list["sns"].lineplot, x, y)
                elif plot_type == "ecdf":   
                        g = plotter_list["sns"].FacetGrid(data, col_wrap=col_wrap, col=col, hue=hue, aspect=1).map(plotter_list["sns"].ecdfplot, x)
                elif plot_type == "lmplot":
                        g = plotter_list["sns"].lmplot(data=data, x=x, y=y, hue=hue, col=col, col_wrap=col_wrap)

                if x_label: g.set_xlabels(x_label)
                if y_label: g.set_ylabels(y_label)
                g.add_legend()
                g.set_titles(col_template="{col_name}")
                if suptitle is not None:
                        g.fig.subplots_adjust(top=top)
                        g.fig.suptitle(suptitle,fontsize=20)

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

        def modify_by_cols(file, ncols, mod):
                mod_file = []
                mod_file.append(plotter.hash_vars[file][0])
                for idx, row in enumerate(plotter.hash_vars[file][1:]):
                        mod_row = row
                        for col in ncols:
                                mod_row[col] = mod(mod_row[col])     
                        mod_file.append(mod_row)
                return mod_file

        for table in plotter.hash_vars.keys():
                if table == "parsed_non_integrated_rank_summary" or table == "parsed_integrated_rank_summary":
                        parse_table(table, include_header=True)
                else:
                        parse_table(table)
        if plotter.hash_vars.get('parsed_non_integrated_rank_summary') is not None:
                order_columns('parsed_non_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_integrated_rank_summary') is not None:
                order_columns('parsed_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_non_integrated_rank_pos_cov') is not None:
                order_columns('parsed_non_integrated_rank_pos_cov',0)
                plotter.hash_vars["parsed_non_integrated_rank_pos_cov"] = modify_by_cols("parsed_non_integrated_rank_pos_cov", [3], lambda x: float(x)/number_of_positives * 100)
                plotter.hash_vars['parsed_non_integrated_rank_pos_cov'] = parse_heatmap_from_flat(plotter.hash_vars['parsed_non_integrated_rank_pos_cov'][1:],1,2,3)

        if plotter.hash_vars.get('parsed_integrated_rank_pos_cov') is not None:
                order_columns('parsed_integrated_rank_pos_cov',0)
                plotter.hash_vars["parsed_integrated_rank_pos_cov"] = modify_by_cols("parsed_integrated_rank_pos_cov", [3], lambda x: float(x)/number_of_positives * 100)
                plotter.hash_vars['parsed_integrated_rank_pos_cov'] = parse_heatmap_from_flat(plotter.hash_vars['parsed_integrated_rank_pos_cov'][1:],1,2,3)

        if plotter.hash_vars.get('parsed_annotation_grade_metrics') is not None:
                order_columns('parsed_annotation_grade_metrics',0)

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
${make_title("figure", "seed_wflow", """Workflow of the benchmarking process. The reference for gene backup controls have been 
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
${make_title("figure", "coverage_bars", """Coverage obtained in each individual (A)
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
${make_title("figure", "rank_boxplot", f"""Rank distributions in each individual (A)
 or integrated (B) eGSM. In both plots, y axis ({italic("Normalized ranks")}) represent the rank is normalized on 0-1 range.""")}

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
<% text.append(make_title("figure", "agg_tops", f"""Top 5,10,20,100 on different individual (A) and integrated (B) eGSM."""))%>
${collapsable_data("Tops absolute values", None, "tops_absolute", "\n".join(text))}

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
${make_title("figure", "agg_tops", f"""Top 5,10,20,100 on different individual (A) and integrated (B) eGSM.""")}
<% txt="Performance curves" %>
${plotter.create_title(txt, id='perf_curves', hlevel=2, indexable=True, clickable=False)}

<% txt="CDF" %>
${plotter.create_title(txt, id='cdf_curves', hlevel=3, indexable=True, clickable=False)}

<a href="https://academic.oup.com/bib/article/23/2/bbac019/6521702#330302198">Xiao Yuan et al. Evaluation of phenotype-driven gene prioritization methods for Mendelian diseases, Briefings in Bioinformatics, Volume 23, Issue 2, March 2022, bbac019 </a>
<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="non_integrated_rank_cdf", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="ecdf",data=data, 
                                        plotter_list=plotter_list, x="rank", col="annot", 
                                        hue="Embedding", col_wrap=3, 
                                        suptitle="A", x_label="Normalized Rank", y_label="TPR", top=0.9))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="integrated_rank_cdf", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="ecdf",data=data, plotter_list=plotter_list, x="rank", 
                                        col="integration", hue="Embedding", col_wrap=2, suptitle="B", x_label="Normalized Rank", y_label="TPR", top=0.8))}
        % endif
</div>
${make_title("figure", "cdf_curve", f"""CDF curves by each individual (A)
 or integrated (B) eGSM. In both plots, y axis represent 
 the true positive rate ({italic("TPR")}) and x axis ({italic("Normalized Rank")}) the rank normalized from 0 to 1.""")}

<% txt="ROC" %>
${plotter.create_title(txt, id='roc_curves', hlevel=3, indexable=True, clickable=False)}

<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_measures") is not None: 

                 ${ plotter.static_plot_main( id="non_integrated_rank_measures", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="lineplot", data=data,
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='annot', 
                                        hue='Embedding', col_wrap=3, suptitle="A", 
                                        top=0.9, x_label="FPR", y_label="TPR"))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_measures") is not None: 
                 ${ plotter.static_plot_main( id="integrated_rank_measures", header=True, row_names=False, smp_attr=[0,1,2,3], fields =[4,5,6], 
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="lineplot",data=data, 
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='integration', 
                                        hue='Embedding', col_wrap=2, suptitle="B", 
                                        top=0.8, labels = 'Embedding', x_label="FPR", y_label="TPR"))}
        % endif
</div>
${make_title("figure", "roc_curve", f"""ROC in each individual (A) or integrated (B) eGSM.""")}

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
                        'setMinX': 0,
                        'setMaxX': 1,
                        "smpTextRotate": 45,
                        "segregateSamplesBy": "Integration"
                        })) %>
% endif
<% txt.append(make_title("figure", "roc_ic", f"""AUROC Confidence Interval (CI) in each individual (A) or integrated (B) eGSM. CI was obtained by a 1000 iteration bootstrap.""")) %>
${collapsable_data("AUROC Confidence Interval", None, "auroc_ci", "\n".join(txt))}