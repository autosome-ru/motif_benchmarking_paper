var count_elements = function(collection) {
	counts = {}
	collection.forEach(el => counts[el] = (counts[el] || 0) + 1);
	return counts;
};
var uniq = function(collection) {
	return [...new Set(collection)];
};

var url_string = window.location.href;
var url = new URL(url_string);
var tsv_fn = url.searchParams.get("tsv");
data = d3.tsv(tsv_fn, d3.autoType);

tfclass_id = function(str){
	var from = str.indexOf('{');
	var to = str.indexOf('}');
	if (from >= 0 && to >= 0) {
		return str.substr(from + 1, to - from - 1 ).split('.').map(x => Number(x));
	} else {
		return [Infinity]; // motifs without TFClass should go in the end
	}
}
tfclass_comparator = function(a, b) {
	var a_fam = tfclass_id(a);
	var b_fam = tfclass_id(b);
	if (a_fam < b_fam) {
	    return -1;
	} else if (a_fam > b_fam) {
	    return 1;
	}
	return 0;
}

prepare_families_graph = function(data) {
	var fam2fam = data.map(function(row){
		// return {from: row['tfs_of_best_motif'], to: row['experiment_TF_family']};
		return {from: row['best_motifs_family'], to: row['experiment_TF_family']};
	});
	fam_counts = count_elements(fam2fam.map(o => [o.from, o.to]).flat());
	fams = Object.keys(fam_counts);
	fams_sorted = Object.keys(fam_counts).sort(function(a,b) {return fam_counts[b] - fam_counts[a]; });
	src_nodes = uniq( fam2fam.map(fam => fam.from + ':from') ).sort(tfclass_comparator).map(name => Object({name: name}));
	dst_nodes = uniq( fam2fam.map(fam => fam.to + ':to') ).sort(tfclass_comparator).map(name => Object({name: name}));
	nodes = src_nodes.concat(dst_nodes);
	links = [];
	for (let src_fam in fams) {
		for (let dst_fam in fams) {
			var count = fam2fam.filter(row => (row.from == fams[src_fam]) && (row.to == fams[dst_fam])).length;
			if (count > 0) {
				link = {
					source: fams[src_fam] + ':from',
					target: fams[dst_fam] + ':to',
					value: count,
					names: [fams[src_fam], fams[dst_fam]],
				};
				links.push(link);
			}
		}
	}

	graph = {
		nodes: nodes,
	 	links: links
	};
	return graph;
}

draw_sankey = function(graph) {
	var width = 1200, height = 650;
	var color = d3.scaleOrdinal(
		// ["other", "unknown"],
		// ["#AAAAAA","#4D4D4D","#F15854","#FAA43A","#e5d00d","#60BD68","#5DA5DA","#F17CB0","#975597","#B2912F",]
		[
			'C2H2 zinc finger factors{2.3}', // several zinc finger families are here
			'Fork head / winged helix factors{3.3}', // FOX-s are here
			'Tryptophan cluster factors{3.5}', // ETS-s are here

			'Three-zinc finger Kruppel-related factors{2.3.1}',
			'More than 3 adjacent zinc finger factors{2.3.3}',
			'Forkhead box (FOX) factors{3.3.1}', // FOX-s are here
			'Ets-related factors{3.5.2}', // ETS-s are here
		],
		[
			'#F15854',
			'#60BD68',
			'#5DA5DA',

			'#faa43a',
			'#F15854',
			'#60BD68',
			'#5DA5DA',
		]
	).unknown("#DDDDDD");

	sankey = d3.sankey()
	    .nodeSort(null)
	    // .linkSort(null)
	    .nodeId(d => d.name)
	    .nodeWidth(4)
	    .iterations(1000)
	    // .nodePadding(20)
	    .extent([[350, 10], [width - 250, height - 10]]);
	graph = sankey(graph);

	var svg = d3.select('body').append('svg').attr('id', 'svg')
		.attr('width', width)
		.attr('height', height)
		.style('width', '100%')
		.style('height', 'auto');

	svg.append("g")
	.selectAll("rect")
	.data(graph.nodes)
	.join("rect")
	  .attr("x", d => d.x0)
	  .attr("y", d => d.y0)
	  .attr("height", d => d.y1 - d.y0)
	  .attr("width", d => d.x1 - d.x0)
	.append("title")
	  .text(d => `${d.name}\n${d.value.toLocaleString()}`);

	svg.append("g")
	  .attr("fill", "none")
	.selectAll("g")
	.data(graph.links)
	.join("path")
	  .attr("d", d3.sankeyLinkHorizontal())
	  .attr("stroke", d => color(d.names[0].replace('\u00FC', 'u') )) // Krüppel --> Kruppel
	  // .attr("stroke", d => color(d.names[1].replace('\u00FC', 'u') )) // Krüppel --> Kruppel
	  .attr("stroke-width", d => d.width)
	  .style("mix-blend-mode", "multiply")
	.append("title")
	  .text(d => `${d.names.join(" → ")}\n${d.value.toLocaleString()}`);

	svg.append("g")
	  .style("font", "10px sans-serif")
	.selectAll("text")
	.data(graph.nodes)
	.join("text")
		.attr("font-family", "Lato")
	  .attr("x", d => d.x0 < width / 2 ? d.x1 - 10 : d.x0 + 10)
	  .attr("y", d => (d.y1 + d.y0) / 2)
	  .attr("dy", "0.35em")
	  .attr("text-anchor", d => d.x0 < width / 2 ? "end" : "start")
	  .text(d => d.name.replace(/:from$/, '').replace(/:to$/, ''))
	.append("tspan")
	  .attr("fill-opacity", 0.7)
	  .text(d => ` ${d.value.toLocaleString()}`);
}

// from https://stackoverflow.com/a/23218877/10712892
add_export_button = function() {
	var svg_node = document.getElementById("svg");

	var serializer = new XMLSerializer();
	var source = serializer.serializeToString(svg_node);

	//add name spaces.
	if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){
	    source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
	}
	if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){
	    source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
	}

	//add xml declaration
	source = '<?xml version="1.0" standalone="no"?>\r\n' + source;

	//convert svg source to URI data scheme.
	var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);

	//set url value to a element's href attribute.
	//you can download svg file by right click menu.
	var link_elem = document.createElement('a');
	link_elem.href = url;
	var link_content = document.createTextNode("Download SVG");
	link_elem.appendChild(link_content);
	document.body.appendChild(link_elem);
}

data.then(function(data){
	graph = prepare_families_graph(data);
	draw_sankey(graph);
	add_export_button();
});
