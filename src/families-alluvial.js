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
	for(var i = 0; i < a_fam.length; ++i) {
		if (a_fam[i] < b_fam[i]) {
			return -1;
		} else if (a_fam[i] > b_fam[i]) {
			return 1;
		}
	}
	return 0;
}

lexicographic_comparator = function(a, b) {
	var a_upper = a.toUpperCase();
	var b_upper = b.toUpperCase();
	if (a_upper < b_upper) {
		return -1;
	} else if (a_upper > b_upper) {
		return 1;
	} else {
		return 0;
	}
}


prepare_families_graph = function(data, key_from, key_to) {
	var fam2fam = data.map(function(row){
		// return {from: row.motif_family + ' - ' + row[key_from], to: row[key_to] + ' - ' + row.experiment_family, num_experiments: row['num_experiments'], ...row};
		return {from: row[key_from], to: row[key_to], num_experiments: row['num_experiments'], ...row};
	});
	// fam2fam = fam2fam.filter(row => row.num_experiments >= 2)
	fam2fam = fam2fam.filter(row => (row.motif_family != 'ambiguous') && (row.experiment_family != 'ambiguous'))
	fam2fam = fam2fam.filter(row => (row.motif_family != 'unknown') && (row.experiment_family != 'unknown'))
	// fam2fam = fam2fam.filter(row => [2, 3] <= tfclass_id(row.motif_family) && tfclass_id(row.motif_family) < [2, 4]);
	// fam2fam = fam2fam.filter(row => [3, 3] <= tfclass_id(row.motif_family) && tfclass_id(row.motif_family) < [3, 4]);
	// fam2fam = fam2fam.filter(row => [3, 5] <= tfclass_id(row.motif_family) && tfclass_id(row.motif_family) < [3, 6]);
	fam_counts = count_elements(fam2fam.map(o => [o.from, o.to]).flat());
	fams = Object.keys(fam_counts);
	var tf_classes = uniq(fams.map(fam => 'separator{' + tfclass_id(fam).slice(0,2).join('.') + '.99}'));
	src_nodes = uniq( fam2fam.map(fam => fam.from).concat(tf_classes).map(fam_name => fam_name + ':from') ).sort(lexicographic_comparator).map(name => Object({name: name}));
	// nodes.push({name: tf_class + ':from'});
	dst_nodes = uniq( fam2fam.map(fam => fam.to).concat(tf_classes).map(fam_name => fam_name + ':to')  ).sort(lexicographic_comparator).map(name => Object({name: name}));
	// nodes.push({name: tf_class + ':to'});
	nodes = src_nodes.concat(dst_nodes);
	links = [];
	// for (let idx in tf_classes) {
	// 	let tf_class = tf_classes[idx];
	// 	link = {
	// 		source: tf_class + ':from',
	// 		target: tf_class + ':to',
	// 		value: 5,
	// 		names: [tf_class, tf_class],
	// 	};
	// 	links.push(link);
	// }
	console.log(fam_counts);
	for (let src_fam in fams) {
		for (let dst_fam in fams) {
			var count = fam2fam.filter(row => (row.from == fams[src_fam]) && (row.to == fams[dst_fam])).length;
			if (count > 0) {
				link = {
					source: fams[src_fam] + ':from',
					target: fams[dst_fam] + ':to',
					// value: count, // for best_motifs
					value: ((count >= 100) ? count : 100) ** 0.5, // for good_motifs
					names: [fams[src_fam], fams[dst_fam]],
				};
				links.push(link);
			}
		}
	}
	console.log(links);

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
			'C2H2 zinc finger factors{2.3}', 'C2H2 ZF', // several zinc finger families are here
			'Fork head / winged helix factors{3.3}', 'Forkhead', // FOX-s are here
			'Tryptophan cluster factors{3.5}', 'Ets', /*'IRF', 'Myb/SANT',*/ // ETS-s are here

			'Three-zinc finger Kruppel-related factors{2.3.1}', ...tfs_by_family['kruppel'],
			'Factors with multiple dispersed zinc fingers{2.3.4}', ...tfs_by_family['dispersed_zinc'],
			'Forkhead box (FOX) factors{3.3.1}', ...tfs_by_family['fox'], // FOX-s are here
			'Ets-related factors{3.5.2}', ...tfs_by_family['ets'], // ETS-s are here
		],
		[
			'#F15854', '#F15854',
			'#60BD68', '#60BD68',
			'#5DA5DA', '#5DA5DA', //'#5DA5DA', '#5DA5DA',

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
	.data(graph.nodes.filter(node => !node.name.startsWith('separator')))
	.join("rect")
	  .attr("x", d => d.x0)
	  .attr("y", d => d.y0)
	  .attr("height", d => d.y1 - d.y0)
	  .attr("width", d => d.x1 - d.x0)
	.append("title")
		.attr("font-family", "Lato")
	  .text(d => `${d.name}\n${d.value.toLocaleString()}`);

	svg.append("g")
	  .attr("fill", "none")
	.selectAll("g")
	.data(graph.links.filter(link => !link.names[0].startsWith('separator')))
	.join("path")
	  .attr("d", d3.sankeyLinkHorizontal())
	  // .attr("stroke", d => color(d.names[color_idx].replace('\u00FC', 'u').split(' - ')[color_idx] )) // Krüppel --> Kruppel
	  .attr("stroke", d => color(d.names[color_idx].replace('\u00FC', 'u') )) // Krüppel --> Kruppel
	  .attr("stroke-width", d => d.width)
	  .style("mix-blend-mode", "multiply")
	.append("title")
	  .text(d => `${d.names.join(" \u2192 ")}\n${d.value.toLocaleString()}`);

	svg.append("g")
	  .style("font", "10px sans-serif")
	.selectAll("text")
	.data(graph.nodes.filter(node => !node.name.startsWith('separator')))
	.join("text")
		.attr("font-family", "Lato")
	  .attr("x", d => d.x0 < width / 2 ? d.x1 - 10 : d.x0 + 10)
	  .attr("y", d => (d.y1 + d.y0) / 2)
	  .attr("dy", "0.35em")
	  .attr("text-anchor", d => d.x0 < width / 2 ? "end" : "start")
	  .text(d => d.name.replace(/:from$/, '').replace(/:to$/, ''))
	// .append("tspan")
	//   .attr("fill-opacity", 0.7)
	//   .text(d => ` ${d.value.toLocaleString()}`);
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
