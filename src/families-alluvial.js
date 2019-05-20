var count_elements = function(collection) {
	counts = {}
	collection.forEach(el => counts[el] = (counts[el] || 0) + 1);
	return counts;
};
var uniq = function(collection) {
	return [...new Set(collection)];
};

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

prepare_families_graph = function(data, key_from, key_to) {
	var fam2fam = data.map(function(row){
		return {from: row[key_from], to: row[key_to]};
	});
	fam_counts = count_elements(fam2fam.map(o => [o.from, o.to]).flat());
	fams = Object.keys(fam_counts);
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

draw_sankey = function(graph, colorize_by, width, height) {
	var tfs_by_family = {
		kruppel: ['EGR1', 'EGR2', 'EGR3', 'EGR4', 'KLF12', 'KLF13', 'KLF14', 'KLF15', 'KLF16', 'KLF1', 'KLF3', 'KLF4', 'KLF5', 'KLF6', 'KLF8', 'KLF9', 'SP1', 'SP2', 'SP3', 'SP4', 'SP8'],
		dispersed_zinc: ['BCL11A', 'E4F1', 'HIC1', 'HIC2', 'HINFP', 'HIVEP1', 'HIVEP2', 'IKZF1', 'INSM1', 'MAZ', 'MECOM', 'PATZ1', 'PRDM4', 'REST', 'RREB1', 'SALL4', 'VEZF1', 'ZBTB17', 'ZBTB4', 'ZNF134', 'ZNF219', 'ZNF335', 'ZNF341', 'ZNF382', 'ZNF418', 'ZNF423', 'ZNF467', 'ZNF770', 'ZNF784', 'ZNF8'],
		fox: ['FOXA1', 'FOXA2', 'FOXA3', 'FOXB1', 'FOXC1', 'FOXC2', 'FOXD1', 'FOXD2', 'FOXD3', 'FOXF1', 'FOXF2', 'FOXG1', 'FOXH1', 'FOXI1', 'FOXJ2', 'FOXJ3', 'FOXK1', 'FOXK2', 'FOXL1', 'FOXM1', 'FOXO1', 'FOXO3', 'FOXO4', 'FOXO6', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXQ1'],
		ets: ['EHF', 'ELF1', 'ELF2', 'ELF3', 'ELF4', 'ELF5', 'ELK1', 'ELK3', 'ELK4', 'ERF', 'ERG', 'ETS1', 'ETS2', 'ETV1', 'ETV2', 'ETV3', 'ETV4', 'ETV5', 'ETV6', 'ETV7', 'FEV', 'FLI1', 'GABPA', 'SPDEF', 'SPI1', 'SPIB', 'SPIC'],
	};

	var color = d3.scaleOrdinal(
		// ["other", "unknown"],
		// ["#AAAAAA","#4D4D4D","#F15854","#FAA43A","#e5d00d","#60BD68","#5DA5DA","#F17CB0","#975597","#B2912F",]
		[
			'C2H2 zinc finger factors{2.3}', // several zinc finger families are here
			'Fork head / winged helix factors{3.3}', // FOX-s are here
			'Tryptophan cluster factors{3.5}', // ETS-s are here

			'Three-zinc finger Kruppel-related factors{2.3.1}', ...tfs_by_family['kruppel'],
			'Factors with multiple dispersed zinc fingers{2.3.4}', ...tfs_by_family['dispersed_zinc'],
			'Forkhead box (FOX) factors{3.3.1}', ...tfs_by_family['fox'], // FOX-s are here
			'Ets-related factors{3.5.2}', ...tfs_by_family['ets'], // ETS-s are here
		],
		[
			'#F15854',
			'#60BD68',
			'#5DA5DA',

			...Array(1 + tfs_by_family['kruppel'].length).fill('#faa43a'),
			...Array(1 + tfs_by_family['dispersed_zinc'].length).fill('#F10854'),
			...Array(1 + tfs_by_family['fox'].length).fill('#60BD68'),
			...Array(1 + tfs_by_family['ets'].length).fill('#5DA5DA'),
		]
	).unknown("#DDDDDD");

	var color_idx = 0;
	if (colorize_by == 'from') {
		color_idx = 0;
	} else if (colorize_by == 'to') {
		color_idx = 1;
	}

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
	  .attr("stroke", d => color(d.names[color_idx].replace('\u00FC', 'u') )) // Krüppel --> Kruppel
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
