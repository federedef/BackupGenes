<%
        import os.path
        import warnings
        import pandas as pd
        import re


        warnings.simplefilter(action='ignore', category=FutureWarning)

        # Text
        #######

        def italic(txt):
                return f"<i>{txt}</i>"

        def collapsable_data(click_title, click_id, txt):
                collapsable_txt = f"""
                {plotter.create_title(click_title, id=None, indexable=False, clickable=True, t_id=click_id)}\n
                <div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
                        {plotter.create_collapsable_container(click_id, txt)}
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
                        plotter.hash_vars[name].insert(0, tab_header)
                else:
                        plotter.hash_vars[name] = parse_data(plotter.hash_vars[name])

        for table in plotter.hash_vars.keys():
                parse_table(table)
        if plotter.hash_vars.get('parsed_non_integrated_rank_summary') is not None:
                order_columns('parsed_non_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_integrated_rank_summary') is not None:
                order_columns('parsed_integrated_rank_summary',0)

        if plotter.hash_vars.get('parsed_non_integrated_rank_pos_cov') is not None:
                order_columns('parsed_non_integrated_rank_pos_cov',0)

        if plotter.hash_vars.get('parsed_integrated_rank_pos_cov') is not None:
                order_columns('parsed_integrated_rank_pos_cov',0)

        if plotter.hash_vars.get('parsed_annotation_grade_metrics') is not None:
                order_columns('parsed_annotation_grade_metrics',0)

        img_path="/mnt/scratch/users/bio_267_uma/federogc/executions/BackupGenes/report/img/"

        
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

        def get_medianrank_size(var_name, groupby = ['annot_kernel','annot','kernel'], value = 'absolute_ranking', tops = [5,10,20]):
                df = pd.DataFrame(plotter.hash_vars[var_name][1:], columns = plotter.hash_vars[var_name][0])
                df_final = df.groupby(groupby)[value].size().reset_index()
                for top in tops:
                        df_top = df.groupby(groupby)[value].apply(lambda x: sum(pd.to_numeric(x)<=top)).reset_index()
                        df_top.rename(columns={value: f"Top {top}"}, inplace=True)
                        df_final = pd.concat([df_final, df_top[[f"Top {top}"]]], axis=1)
                print(df_final.columns)
                #col_names = plotter.hash_vars[var_name][0]
                #col_names.append("size")
                return [df_final.columns.tolist()] + df_final.values.tolist()


        plotter.hash_vars["non_integrated_tops"] = get_medianrank_size("non_integrated_rank_cdf", tops=[5,10,20,100])
        plotter.hash_vars["integrated_tops"] = get_medianrank_size("integrated_rank_cdf", ["integration_kernel","integration","kernel"], tops = [5,10,20,100])
        
%>

<% plotter.set_header() %>

<% txt="From eGSM to Backup Genes." %>
${plotter.create_title(txt, id='main_backup_gene', hlevel=1, indexable=True, clickable=False)}


<% txt="Workflow for Backup Gene Benchmarking" %>
${plotter.create_title(txt, id='main_backup_gene', hlevel=2, indexable=True, clickable=False)}

<%
        graph=f"""
        ---
        title: Backup Benchmarking Flux
        config:
         theme: dark
         themeVariables:
          lineColor: "#717171"
        
        ---
        graph LR
         SI[Manually curated from literature]
         SII[Double knock-out \\n screening]
         Papi[<span style="color:#280054">Big Papi \\n screening</span>]
         Digenic[<span style="color:#280054">Digenic \\n screening</span>]
         P[<span style="color:#023020">Positive</span>]
         N[<span style="color:#8B0000">Negative</span>]
         CDF[<span style="color:#000000">CDF</span>]
         ROC[<span style="color:#000000">ROC with Bootstrap \\n 1000 iterations</span>]
         Coverage[<span style="color:#000000">Coverage</span>]
         subgraph C[Backup Control \\n Dataset]
          P
          N
         end
         SI--> P
         SII-->Papi
         SII-->Digenic
         Papi -- Positive or negative interaction \\n based on p-value = .05 cutoff \\n among different cell lines --> C
         Digenic -- Positive or negative interaction \\n based on p-value = .05 cutoff \\n among different cell lines --> C
         P--> CDF
         P--> Coverage
         P & N --> ROC
         style P fill:#84D677
         style N fill:#FF503E
         style CDF fill:#A0A0A0
         style ROC fill:#A0A0A0
         style Coverage fill:#A0A0A0
         style Papi fill:#BF7AE7
         style Digenic fill:#BF7AE7
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
                        ${plotter.barplot(id='parsed_non_integrated_rank_pos_cov', responsive= False, header=True,
                         fields = [1,3],
                         x_label = 'Number of control candidate \n genes present',
                         height = '400px', width= '400px',
                         var_attr = [1,2],
                         title = "(A) Individual eGSM",
                         config = {
                                'showLegend' : True,
                                'graphOrientation' : 'horizontal',
                                'colorBy' : 'Kernel',
                                'setMinX': 0,
                                "titleFontStyle": "italic",
                                "titleScaleFontFactor": 0.3
                                })}
        % endif
        </div>
        <div style="margin-left: 10px;"> 
                % if plotter.hash_vars.get('parsed_integrated_rank_pos_cov') is not None: 
                        ${plotter.barplot(id= "parsed_integrated_rank_pos_cov", fields= [1,3] , header= True, responsive= False,
                                height= '400px', width= '400px', x_label= 'Number of control candidate \n genes present' , var_attr= [1,2],
                                title = "(B) Integrated eGSM",
                                config = {
                                        'showLegend' : True,
                                        'graphOrientation' : 'horizontal',
                                        'colorBy' : 'Kernel',
                                        'setMinX': 0,
                                        "titleFontStyle": "italic",
                                        "titleScaleFontFactor": 0.3
                                        })}
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
                                ${plotter.boxplot(id= 'non_integrated_rank_cdf', header= True, row_names= False, default= False, fields= [5],  var_attr= [0,1,2], group = "kernel",
                                   title= "(A) Individual eGSM",
                                        x_label= "Normalized rank",
                                        config= {
                                                "graphOrientation": "vertical",
                                                "colorBy" : "kernel",
                                                "groupingFactors" :
                                                ["kernel"],
                                                "titleFontStyle": "italic",
                                                "titleScaleFontFactor": 0.3,
                                                "segregateSamplesBy": "annot"})}
                        % endif
        </div>
        <div style="margin-left: 10px;">
                        % if plotter.hash_vars.get('integrated_rank_cdf') is not None: 
                                ${plotter.boxplot(id= 'integrated_rank_cdf', header= True, row_names= False, default= False, fields = [5], var_attr= [0,1,2], group= "kernel", 
                                        title= "(B) Integrated eGSM",
                                        xlabel= "Normalized rank",
                                        config= {
                                                "graphOrientation": "vertical",
                                                "colorBy" : "kernel",
                                                "xAxisTitle": "Normalized rank",
                                                "groupingFactors" :
                                                ["kernel"],
                                                "titleFontStyle": "italic",
                                                "titleScaleFontFactor": 0.3,
                                                "segregateSamplesBy": "integration"})}
                        % endif
        </div>
</div>
${make_title("figure", "rank_boxplot", f"""Rank distributions in each individual (A)
 or integrated (B) eGSM. In both plots, y axis ({italic("Normalized ranks")}) represent the rank is normalized on 0-1 range.""")}

<% txt="Tops" %>
${plotter.create_title(txt, id='tops', hlevel=3, indexable=True, clickable=False)}

<div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
        ${plotter.barplot(id="non_integrated_tops", header=True, fields=[2,4,5,6,7], var_attr=[1,2], 
                title= "(A) Individual eGSM", 
                x_label="Number of top-ranked Backups",
                config={
                  'graphOrientation' : 'vertical',
                  'segregateSamplesBy': "annot",
                  "titleFontStyle": "italic",
                 'setMaxX': 80,
                 "smpLabelRotate": 45
                })}
        ${plotter.barplot(id="integrated_tops", header=True, fields=[2,4,5,6,7], var_attr=[1,2], 
                title= "(B) Integrated eGSM",
                x_label="Number of top-ranked Backups", 
                config={
                  'graphOrientation' : 'vertical',
                  'segregateSamplesBy': "integration",
                  "titleFontStyle": "italic",
                  "titleScaleFontFactor": 0.2,
                 "smpLabelRotate": 45,
                 "smpLabelScaleFontFactor": 0.3,
                 'setMaxX': 80,
                })}
</div>
${make_title("figure", "agg_tops", f"""Top 5,10,20,100 on different individual (A) and integrated (B) eGSM.""")}
<% txt="Performance curves" %>
${plotter.create_title(txt, id='perf_curves', hlevel=2, indexable=True, clickable=False)}

<% txt="CDF" %>
${plotter.create_title(txt, id='cdf_curves', hlevel=3, indexable=True, clickable=False)}
<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="non_integrated_rank_cdf", header=True, row_names=False, var_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="ecdf",data=data, 
                                        plotter_list=plotter_list, x="rank", col="annot", 
                                        hue="kernel", col_wrap=4, 
                                        suptitle="A", x_label="Normalized Rank", y_label="TPR", top=0.9))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_cdf") is not None: 
                ${ plotter.static_plot_main( id="integrated_rank_cdf", header=True, row_names=False, var_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="ecdf",data=data, plotter_list=plotter_list, x="rank", 
                                        col="integration", hue="kernel", col_wrap=2, suptitle="B", x_label="Normalized Rank", y_label="TPR", top=0.8))}
        % endif
</div>
${make_title("figure", "cdf_curve", f"""CDF curves by each individual (A)
 or integrated (B) eGSM. In both plots, y axis represent 
 the true positive rate ({italic("TPR")}) and x axis ({italic("Normalized Rank")}) the rank normalized from 0 to 1.""")}

<% txt="ROC" %>
${plotter.create_title(txt, id='roc_curves', hlevel=3, indexable=True, clickable=False)}

<div style="overflow: hidden; text-align:center">
        % if plotter.hash_vars.get("non_integrated_rank_measures") is not None: 
                 ${ plotter.static_plot_main( id="non_integrated_rank_measures", header=True, row_names=False, var_attr=[0,1,2,3], fields =[4,5,6],
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="lineplot", data=data,
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='annot', 
                                        hue='kernel', col_wrap=4, suptitle="A", 
                                        top=0.9, x_label="FPR", y_label="TPR"))}
        % endif
        % if plotter.hash_vars.get("integrated_rank_measures") is not None: 
                 ${ plotter.static_plot_main( id="integrated_rank_measures", header=True, row_names=False, var_attr=[0,1,2,3], fields =[4,5,6], 
                                plotting_function= lambda data, plotter_list: plot_with_facet(plot_type="lineplot",data=data, 
                                        plotter_list=plotter_list, x='fpr', y='tpr', col='integration', 
                                        hue='kernel', col_wrap=2, suptitle="B", 
                                        top=0.8, labels = 'kernel', x_label="FPR", y_label="TPR"))}
        % endif
</div>
${make_title("figure", "roc_curve", f"""ROC in each individual (A) or integrated (B) eGSM.""")}

<%txt=[]%>
% if plotter.hash_vars.get("parsed_non_integrated_rank_summary") is not None: 
        <% txt.append(plotter.line(id= "parsed_non_integrated_rank_summary", fields= [0, 7, 13, 8], header= True, row_names= True,
                responsive= False,
                height= '400px', width= '400px', x_label= 'AUC',
                title= "(A) Individual eGSM",
                config= {
                        'showLegend' : True,
                        'graphOrientation' : 'vertical',
                        "titleFontStyle": "italic",
                        "titleScaleFontFactor": 0.3,
                        'setMinX': 0,
                        'setMaxX': 1,
                        "smpLabelRotate": 45
                        })) %>
% endif
% if plotter.hash_vars.get('parsed_integrated_rank_summary') is not None: 
        <% txt.append(plotter.line(id= "parsed_integrated_rank_summary", fields=  [0, 7, 13, 8], header= True, row_names= True,
                responsive= False,
                height= '400px', width= '400px', x_label= 'AUC',
                title= "(B) Integrated eGSM",
                config= {
                        'showLegend' : True,
                        'graphOrientation' : 'vertical',
                        "titleFontStyle": "italic",
                        "titleScaleFontFactor": 0.3,
                        'setMinX': 0,
                        'setMaxX': 1,
                        "smpLabelRotate": 45
                        })) %>
% endif
<% txt.append(make_title("figure", "roc_ic", f"""AUROC Confidence Interval (CI) in each individual (A) or integrated (B) eGSM. CI was obtained by a 1000 iteration bootstrap.""")) %>
${collapsable_data("AUROC Confidence Interval", "auroc_ci", "\n".join(txt))}