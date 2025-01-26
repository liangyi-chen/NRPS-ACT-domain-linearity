# -*- coding: utf-8 -*-
"""
Created on Sun May 28 17:01:37 2023

@author: rjche
"""

import sys
import os.path
import glob,os
import re
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplcursors
import plotly.graph_objects as go
from plotly.offline import plot

def check_file_for_nrp(file_path):
    with open(file_path) as file:
        data = json.load(file)
        cluster_data = data.get("cluster")
        if cluster_data:
            biosyn_class = cluster_data.get("biosyn_class")
            if biosyn_class and "NRP" in biosyn_class:
                return True
    return False

def check_file_for_pks(file_path):
    with open(file_path) as file:
        data = json.load(file)
        cluster_data = data.get("cluster")
        if cluster_data:
            biosyn_class = cluster_data.get("biosyn_class")
            if biosyn_class and "Polyketide" in biosyn_class:
                return 'NRPS+PKS'
    return 'NRPS'


def extract_info_from_file(file_path):
    with open(file_path) as file:
        data = json.load(file)
        cluster_data = data.get("cluster")
        if cluster_data:
            compounds = cluster_data.get("compounds")
            if compounds and len(compounds) > 0:
                first_compound = compounds[0]
                compound_name = first_compound.get("compound")
                mol_mass = first_compound.get("mol_mass")
                return compound_name, mol_mass
    return None, None

def count_domains(domain_list):
    a_domain_count = 0
    c_domain_count = 0
    t_domain_count = 0

    for domain in domain_list:
        if 'AMP-binding' in domain or 'A-OX' in domain:
            a_domain_count += 1
        #if 'Condensation' in domain or 'Condensation_DCL' in domain or 'Condensation_Dual' in domain or 'Condensation_LCL' in domain or 'Condensation_Starter' in domain or 'Cglyc' in domain:
        if 'Condensation' in domain or 'Heterocyclization' in domain:
            c_domain_count += 1
        if 'PCP' in domain or 'PP-binding' in domain:
            t_domain_count += 1

    return a_domain_count, c_domain_count, t_domain_count

NRP_file = []
for file in glob.glob(('MIBiG_data/*'+'.json')):
    if check_file_for_nrp(file):
        NRP_file.append(file)


'''
NRP_file = ['BGC0002624.json']
'''
csv_output = []
for i in range(len(NRP_file)):
    csv_row = []
    csv_row.append(NRP_file[i].removesuffix('.json'))
    csv_row.append(check_file_for_pks(NRP_file[i]))
    
    with open(NRP_file[i].replace('.json', '.gbk'),"r") as filegbk:
        gbk_content = filegbk.read()
        gbk_content = gbk_content.replace('\n', '')
    domain_matches = re.findall(r'/aSDomain="(.*?)"', gbk_content)
    '''
    if 'PKS_PP' in domain_matches or 'PKS_KS' in domain_matches or 'PKS_AT' in domain_matches or 'PKS_Docking_Cterm' in domain_matches or 'PKS_KR' in domain_matches or 'PKS_DH' in domain_matches or 'PKS_Docking_Nterm' in domain_matches or 'PKS_ER' in domain_matches or 'PKS_DH2' in domain_matches or 'PKS_DHt' in domain_matches:
        if csv_row[1]=='NRPS':
            csv_row = [sub.replace('NRPS', 'wrongNRPS+PKS') for sub in csv_row]
    '''
    csv_row.append(extract_info_from_file(NRP_file[i])[0])
    csv_row.append(extract_info_from_file(NRP_file[i])[1])
    csv_row.append(count_domains(domain_matches)[0])
    csv_row.append(count_domains(domain_matches)[1])
    csv_row.append(count_domains(domain_matches)[2])
    #csv_row.append(domain_matches)
    csv_output.append(csv_row)
csv_output.append(['BGC000xxxx','NRPS', 'fimsbactin', 574.5, 2, 4, 4])
print(csv_output)



output_directory = 'Output/'
# Create the directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = 'output.csv'

# Write the data to a CSV file
with open(output_directory+output_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(csv_output)


filtered_csv = [row for row in csv_output if row[3] is not None and row[1] == 'NRPS'] #remove the rows where Mw is None and type is not NRPS

A_domain_online = [row for row in filtered_csv if float(row[3])>=(57.05132*float(row[4])-18*(float(row[4])-1)) and float(row[3])<=(186.2099*float(row[4])-18*(float(row[4])-1))]
A_domain_outlier = [row for row in filtered_csv if float(row[3])<(57.05132*float(row[4])-18*(float(row[4])-1)) or float(row[3])>(186.2099*float(row[4])-18*(float(row[4])-1))]
C_domain_online = [row for row in filtered_csv if float(row[3])>=(57.05132*(float(row[5])+1)-18*float(row[5])) and float(row[3])<=(186.2099*(float(row[5])+1)-18*float(row[4]))]
C_domain_outlier = [row for row in filtered_csv if float(row[3])<(57.05132*(float(row[5])+1)-18*float(row[5])) or float(row[3])>(186.2099*(float(row[5])+1)-18*float(row[4]))]
T_domain_online = [row for row in filtered_csv if float(row[3])>=(57.05132*float(row[6])-18*(float(row[6])-1)) and float(row[3])<=(186.2099*float(row[6])-18*(float(row[6])-1))]
T_domain_outlier = [row for row in filtered_csv if float(row[3])<(57.05132*float(row[6])-18*(float(row[6])-1)) or float(row[3])>(186.2099*float(row[6])-18*(float(row[6])-1))]
#print(len(filtered_csv),len(T_domain_online),len(T_domain_outlier))

#Do linear fit
xa_online = np.array([row[3] for row in A_domain_online])
ya_online = np.array([row[4] for row in A_domain_online])
aa,ba = np.polyfit(xa_online,ya_online,1)
xa_offline = np.array([row[3] for row in A_domain_outlier])
ya_offline = np.array([row[4] for row in A_domain_outlier])

xc_online = np.array([row[3] for row in C_domain_online])
yc_online = np.array([row[5] for row in C_domain_online])
ac,bc = np.polyfit(xc_online,yc_online,1)
xc_offline = np.array([row[3] for row in C_domain_outlier])
yc_offline = np.array([row[5] for row in C_domain_outlier])

xt_online = np.array([row[3] for row in T_domain_online])
yt_online = np.array([row[6] for row in T_domain_online])
at,bt = np.polyfit(xt_online,yt_online,1)
xt_offline = np.array([row[3] for row in T_domain_outlier])
yt_offline = np.array([row[6] for row in T_domain_outlier])

'''
# Create a scatter plot
fig1, ax = plt.subplots()
scatter_online = ax.scatter(xa_online, ya_online)

# Add annotations to the data points
annotations_online = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in A_domain_online]
cursor = mplcursors.cursor(scatter_online, hover=True)
cursor.connect("add", lambda sel: sel.annotation.set_text(annotations_online[sel.target.index]))

# Configure the plot
ax.set_xlabel('Mw')
ax.set_ylabel('A domain number')
ax.set_title('A domain')


# Save the figure
plt.savefig("scatter_plot.png")
'''
annotations_a_online = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in A_domain_online]
annotations_a_offline = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in A_domain_outlier]
annotations_c_online = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in C_domain_online]
annotations_c_offline = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in C_domain_outlier]
annotations_t_online = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in T_domain_online]
annotations_t_offline = [f'{row[0]}, {row[1]}, {row[2]}, {row[3]}' for row in T_domain_outlier]

fig1=go.Figure()
fig1.add_trace(go.Scatter(x=xa_online, y=ya_online, mode='markers', text=annotations_a_online,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
color='blue',  # Set the marker color
        size=10  # Set the marker size
)))

fig1.add_trace(go.Scatter(x=xa_offline, y=ya_offline, mode='markers', text=annotations_a_offline,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
        color='red',  # Set the marker color
        size=5  # Set the marker size
)))
fig1.add_trace(go.Scatter(x=xa_online, y=aa*xa_online+ba, mode='lines',hoverinfo='none'))
fig1.update_layout(
    xaxis_title="Molecular weight",
    yaxis_title="A domain count",
    showlegend=False,
    plot_bgcolor='white',  # Set the plot background color
    xaxis=dict(
        linecolor='black',  # Set the x-axis line color
        linewidth=2,  # Set the x-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the x-axis grid color
        gridwidth=0.5,  # Set the x-axis grid line width
    ),
    yaxis=dict(
        linecolor='black',  # Set the y-axis line color
        linewidth=2,  # Set the y-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the y-axis grid color
        gridwidth=0.5,  # Set the y-axis grid line width
    )
)
# Generate the HTML file
html_fig1 = '''
<html>
<head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        .input-box {
            margin-bottom: 10px;
            position: relative;
        }
        .message-box {
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 10px;
        }
        .red-text {
            color: red;
        }
    </style>
</head>
<body>
    <div class="message-box" id="messageBox"></div>
    <div class="input-box">
        Custom molecular weight: <input type="number" id="number1" onchange="updatePlot()">
    </div>
    <div class="input-box">
        Custom A domain count: <input type="number" id="number2" onchange="updatePlot()">
    </div>
    <div id="graph"></div>
    <script type="text/javascript">
        var figure = ''' + fig1.to_json() + ''';
        figure.layout.width = 600;  // Set the width of the plot (in pixels)
        figure.layout.height = 400;  // Set the height of the plot (in pixels)
        Plotly.newPlot('graph', figure.data, figure.layout);

        var traceIndex = 3;  // Index of the interactive point trace

        function updatePlot() {
            var xValue = document.getElementById('number1').value;
            var yValue = document.getElementById('number2').value;

            if (xValue && yValue) {
                if (figure.data.length <= traceIndex) {
                    // Create the interactive point trace if it doesn't exist
                    var trace = {
                        x: [],
                        y: [],
                        mode: 'markers',
                        marker: {
                            color: 'orange',
                            size: 15
                        }
                    };
                    figure.data.push(trace);
                }

                // Update the coordinates of the interactive point
                figure.data[traceIndex].x = [parseFloat(xValue)];
                figure.data[traceIndex].y = [parseFloat(yValue)];

                // Update the message box based on the criteria
                var messageBox = document.getElementById('messageBox');
                var x = parseFloat(xValue);
                var y = parseFloat(yValue);
                var linearCondition = x >= 57.05132 * y - 18 * (y - 1) && x <= 186.2099 * y - 18 * (y - 1);
                if (linearCondition) {
                    messageBox.innerHTML = 'This might be a linear NRPS';
                    messageBox.style.color = 'black';
                } else {
                    messageBox.innerHTML = 'This might <span class="red-text">not</span> be a linear NRPS';
                    messageBox.style.color = 'black';
                }

                Plotly.redraw('graph');
            }
        }
    </script>
</body>
</html>
'''
with open(output_directory+"A_domain.html", "w") as file:
    file.write(html_fig1)

fig2=go.Figure()
fig2.add_trace(go.Scatter(x=xc_online, y=yc_online, mode='markers', text=annotations_c_online,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
color='blue',  # Set the marker color
        size=10  # Set the marker size
)))
fig2.add_trace(go.Scatter(x=xc_offline, y=yc_offline, mode='markers', text=annotations_c_offline,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
        color='red',  # Set the marker color
        size=5  # Set the marker size
)))
fig2.add_trace(go.Scatter(x=xc_online, y=ac*xc_online+bc, mode='lines',hoverinfo='none'))
fig2.update_layout(
    xaxis_title="Molecular weight",
    yaxis_title="C domain count",
    showlegend=False,
    plot_bgcolor='white',  # Set the plot background color
    xaxis=dict(
        linecolor='black',  # Set the x-axis line color
        linewidth=2,  # Set the x-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the x-axis grid color
        gridwidth=0.5,  # Set the x-axis grid line width
    ),
    yaxis=dict(
        linecolor='black',  # Set the y-axis line color
        linewidth=2,  # Set the y-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the y-axis grid color
        gridwidth=0.5,  # Set the y-axis grid line width
    )
)
# Generate the HTML file
html_fig2 = '''
<html>
<head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        .input-box {
            margin-bottom: 10px;
            position: relative;
        }
        .message-box {
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 20px;
        }
        .red-text {
            color: red;
        }
    </style>
</head>
<body>
    <div class="message-box" id="messageBox"></div>
    <div class="input-box">
        Custom molecular weight: <input type="number" id="number1" onchange="updatePlot()">
    </div>
    <div class="input-box">
        Custom C domain count: <input type="number" id="number2" onchange="updatePlot()">
    </div>
    <div id="graph"></div>
    <script type="text/javascript">
        var figure2 = ''' + fig2.to_json() + ''';
        figure2.layout.width = 600;  // Set the width of the plot (in pixels)
        figure2.layout.height = 400;  // Set the height of the plot (in pixels)
        Plotly.newPlot('graph', figure2.data, figure2.layout);

        var traceIndex = 3;  // Index of the interactive point trace

        function updatePlot() {
            var xValue = document.getElementById('number1').value;
            var yValue = document.getElementById('number2').value;

            if (xValue && yValue) {
                if (figure2.data.length <= traceIndex) {
                    // Create the interactive point trace if it doesn't exist
                    var trace = {
                        x: [],
                        y: [],
                        mode: 'markers',
                        marker: {
                            color: 'orange',
                            size: 15
                        }
                    };
                    figure2.data.push(trace);
                }

                // Update the coordinates of the interactive point
                figure2.data[traceIndex].x = [parseFloat(xValue)];
                figure2.data[traceIndex].y = [parseFloat(yValue)];

                // Update the message box based on the criteria
                var messageBox = document.getElementById('messageBox');
                var x = parseFloat(xValue);
                var y = parseFloat(yValue);
                var linearCondition = x >= 57.05132 * y - 18 * (y) && x <= 186.2099 * y - 18 * (y);
                if (linearCondition) {
                    messageBox.innerHTML = 'This might be a linear NRPS';
                    messageBox.style.color = 'black';
                } else {
                    messageBox.innerHTML = 'This might <span class="red-text">not</span> be a linear NRPS';
                    messageBox.style.color = 'black';
                }

                Plotly.redraw('graph');
            }
        }
    </script>
</body>
</html>
'''
with open(output_directory+"C_domain.html", "w") as file:
    file.write(html_fig2)

fig3=go.Figure()
fig3.add_trace(go.Scatter(x=xt_online, y=yt_online, mode='markers', text=annotations_t_online,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
color='blue',  # Set the marker color
        size=10  # Set the marker size
)))
fig3.add_trace(go.Scatter(x=xt_offline, y=yt_offline, mode='markers', text=annotations_t_offline,hovertemplate="(%{x}, %{y})<br>%{text}",marker=dict(
        color='red',  # Set the marker color
        size=5  # Set the marker size
)))
fig3.add_trace(go.Scatter(x=xt_online, y=ac*xt_online+bc, mode='lines',hoverinfo='none'))
fig3.update_layout(
    xaxis_title="Molecular weight",
    yaxis_title="T domain count",
    showlegend=False,
    plot_bgcolor='white',  # Set the plot background color
    xaxis=dict(
        linecolor='black',  # Set the x-axis line color
        linewidth=2,  # Set the x-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the x-axis grid color
        gridwidth=0.5,  # Set the x-axis grid line width
    ),
    yaxis=dict(
        linecolor='black',  # Set the y-axis line color
        linewidth=2,  # Set the y-axis line width
        showgrid=True,
        gridcolor='grey',  # Set the y-axis grid color
        gridwidth=0.5,  # Set the y-axis grid line width
    )
)
# Generate the HTML file
html_fig3 = '''
<html>
<head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        .input-box {
            margin-bottom: 10px;
            position: relative;
        }
        .message-box {
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 20px;
        }
        .red-text {
            color: red;
        }
        /* Add these styles to your existing CSS */
        #graph {
            width: 100%; /* Set the width to 100% of its container */
            height: 0; /* Set an initial height to 0, will be adjusted by padding-bottom */
            padding-bottom: 75%; /* Set the aspect ratio (e.g., 75% for 4:3) */
        }
    </style>
</head>
<body>
    <div class="message-box" id="messageBox"></div>
    <div class="input-box">
        Custom molecular weight: <input type="number" id="number1" onchange="updatePlot()">
    </div>
    <div class="input-box">
        Custom T domain count: <input type="number" id="number2" onchange="updatePlot()">
    </div>
    <div id="graph-container">
        <div id="graph"></div>
    </div>
    <script type="text/javascript">
        var figure3 = ''' + fig3.to_json() + ''';
        figure3.layout.width = 600;  // Set the width of the plot (in pixels)
        figure3.layout.height = 400;  // Set the height of the plot (in pixels)
        Plotly.newPlot('graph', figure3.data, figure3.layout);

        var traceIndex = 3;  // Index of the interactive point trace

        function updatePlot() {
            var xValue = document.getElementById('number1').value;
            var yValue = document.getElementById('number2').value;

            if (xValue && yValue) {
                if (figure3.data.length <= traceIndex) {
                    // Create the interactive point trace if it doesn't exist
                    var trace = {
                        x: [],
                        y: [],
                        mode: 'markers',
                        marker: {
                            color: 'orange',
                            size: 15
                        }
                    };
                    figure3.data.push(trace);
                }

                // Update the coordinates of the interactive point
                figure3.data[traceIndex].x = [parseFloat(xValue)];
                figure3.data[traceIndex].y = [parseFloat(yValue)];

                // Update the message box based on the criteria
                var messageBox = document.getElementById('messageBox');
                var x = parseFloat(xValue);
                var y = parseFloat(yValue);
                var linearCondition = x >= 57.05132 * y - 18 * (y-1) && x <= 186.2099 * y - 18 * (y-1);
                if (linearCondition) {
                    messageBox.innerHTML = 'This might be a linear NRPS';
                    messageBox.style.color = 'black';
                } else {
                    messageBox.innerHTML = 'This might <span class="red-text">not</span> be a linear NRPS';
                    messageBox.style.color = 'black';
                }

                Plotly.redraw('graph');
            }
        }
    </script>
</body>
</html>
'''
with open(output_directory+"T_domain.html", "w") as file:
    file.write(html_fig3)



