import py_exp_calc.exp_calc as pxc
import sys
import re
# Text
#######
def italic(txt):
        return f"<i>{txt}</i>"

def bold(txt):
	return f"<b>{txt}</b>"

def collapsable_data(plotter, click_title, click_id, container_id, txt, indexable=False, hlevel=1):
        collapsable_txt = f"""
        {plotter.create_title(click_title, id=click_id, indexable=indexable, clickable=True, hlevel=hlevel, t_id=container_id)}\n
        <div style="overflow: hidden; display: flex; flex-direction: row; justify-content: center;">
                {plotter.create_collapsable_container(container_id, txt)}
        </div>"""
        return collapsable_txt

def make_title(plotter, type, id, sentence):
        if type == "table":
                key = f"tab:{id}"
                html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {plotter.add_table(key)} </b> {sentence} </p>"
        elif type == "figure":
                key = id
                html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {plotter.add_figure(key)} </b> {sentence} </p>"
        return html_title

def ul(lis):
        txt = '<lu class="body_ul">'
        for li in lis:
                li = li.split(":")
                li[0] = bold(li[0])
                li = ":".join(li)
                txt += f"<li>{li}</li>\n"
        txt += "</lu>"
        return  txt

# PARSING TABLES
################


def parse_heatmap_from_flat(data,nrow,ncol,nvalue,smp_attr,var_attr,scale_factor=1):
        pairs = {}
        sample_attributes = {}
        var_attributes = {}
        for row in data:
                if not pairs.get(row[nrow]):
                        pairs[row[nrow]] = {}
                pairs[row[nrow]][row[ncol]] = row[nvalue]
                if smp_attr: sample_attributes[row[nrow]] = [row[attr] for attr in smp_attr.keys()]
                if var_attr: var_attributes[row[ncol]] = [row[attr] for attr in var_attr.keys()]
        mat, row, col = pxc.to_wmatrix_rectangular(pairs)
        mat = scale_factor*mat
        col_attrs = []
        col_attrs = [i for i in smp_attr.values()] if smp_attr else []
        table = [["-", *col_attrs, *col]]
        # Adding var attributes
        if var_attributes:
            for i, attr_name in enumerate(var_attr.values()):
                table.append([attr_name]+["-"]*(len(table[0])-1))
            for cidx, cid in enumerate(col):
                cidx = cidx + len(col_attrs) + 1
                for attr_idx, attr in enumerate(var_attributes[cid]):
                    attr_idx = attr_idx + 1
                    table[attr_idx][cidx] = attr
        if sample_attributes:
            for idx,elem in enumerate(row): table.append([elem,*sample_attributes.get(elem), *mat[idx,:].tolist()])
        else:
            for idx,elem in enumerate(row): table.append([elem, *mat[idx,:].tolist()])
        return table

def parsed_string(data, blacklist = ["sim"]):
        words = []
        for word in data.split("_"):
                for blackword in blacklist:
                        word = re.sub(blackword,"",word)
                word = word.capitalize()
                words.append(word)
        parsed_data = " ".join(words)
        return parsed_data

def parse_data(table, parse_name, blacklist = ["sim"]):
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

def round_table(table, round_by=2):
        rounded_table = []
        for i,row in enumerate(table):
                rounded_table.append(row)
                for j,data in enumerate(row):
                        if data.replace(".", "").isnumeric():
                                rounded_table[i][j] = str(round(float(data),round_by))
                        else:
                                continue
        return rounded_table
        
def order_columns(plotter, name, column):
        tab_header = plotter.hash_vars[name].pop(0)
        plotter.hash_vars[name].sort(key=lambda x: x[column])
        plotter.hash_vars[name].insert(0, tab_header)

def parse_table(plotter, name, parse_name, blacklist=["sim"], include_header = False):
	if not include_header:
		tab_header = plotter.hash_vars[name].pop(0)
		plotter.hash_vars[name] = parse_data(plotter.hash_vars[name], parse_name)
		plotter.hash_vars[name] = round_table(plotter.hash_vars[name])
		plotter.hash_vars[name].insert(0, tab_header)
	else:
		plotter.hash_vars[name] = parse_data(plotter.hash_vars[name], parse_name)
		plotter.hash_vars[name] = round_table(plotter.hash_vars[name])

def parse_layer_names(table, column, get_names):
        data = []
        for row in table:
                data.append(row)
                data[-1][column] = get_name[data[-1][column]]
        return table

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

def modify_by_cols(plotter, file, ncols, mod):
        mod_file = []
        mod_file.append(plotter.hash_vars[file][0])
        for idx, row in enumerate(plotter.hash_vars[file][1:]):
                mod_row = row
                for col in ncols:
                        mod_row[col] = mod(mod_row[col])     
                mod_file.append(mod_row)
        return mod_file