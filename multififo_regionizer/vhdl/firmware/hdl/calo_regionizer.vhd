library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity calo_regionizer is
    port(
            ap_clk : IN STD_LOGIC;
            ap_rst : IN STD_LOGIC;
            ap_start : IN STD_LOGIC;
            ap_done : OUT STD_LOGIC;
            ap_idle : OUT STD_LOGIC;
            ap_ready : OUT STD_LOGIC;
            newevent : IN STD_LOGIC;
            calo_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_out_0_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_1_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_2_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_3_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_4_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_5_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_6_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_7_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_8_V : OUT STD_LOGIC_VECTOR (63 downto 0);
            calo_out_valid_0 : OUT STD_LOGIC;
            calo_out_valid_1 : OUT STD_LOGIC;
            calo_out_valid_2 : OUT STD_LOGIC;
            calo_out_valid_3 : OUT STD_LOGIC;
            calo_out_valid_4 : OUT STD_LOGIC;
            calo_out_valid_5 : OUT STD_LOGIC;
            calo_out_valid_6 : OUT STD_LOGIC;
            calo_out_valid_7 : OUT STD_LOGIC;
            calo_out_valid_8 : OUT STD_LOGIC;
            newevent_out : OUT STD_LOGIC

    );
end calo_regionizer;

architecture Behavioral of calo_regionizer is
    constant NREGIONS  : natural := NPFREGIONS;
    constant NALLFIFOS : natural := NCALOSECTORS*NCALOFIFOS;
    constant NMERGE2   : natural := NALLFIFOS/2;
    constant NMERGE4   : natural := NALLFIFOS/4;

    signal links_in :       particles(NCALOSECTORS*NCALOFIBERS-1 downto 0) := (others => null_particle);
    signal fifo_in :        particles(NALLFIFOS-1 downto 0);
    signal fifo_in_write :  std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_in_roll  :  std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');

    signal fifo_out :         particles(NALLFIFOS-1 downto 0);
    signal fifo_out_valid :   std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_out_full:     std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');
    signal fifo_out_roll:     std_logic_vector(NALLFIFOS-1 downto 0) := (others => '0');

    signal merged2_out :        particles(NMERGE2-1 downto 0);
    signal merged2_out_valid :  std_logic_vector(NMERGE2-1 downto 0) := (others => '0');
    signal merged2_out_roll:    std_logic_vector(NMERGE2-1 downto 0) := (others => '0');
    signal merged2_out_full:    std_logic_vector(NMERGE2-1 downto 0) := (others => '0');

    signal merged4_out :        particles(NMERGE4-1 downto 0);
    signal merged4_out_valid :  std_logic_vector(NMERGE4-1 downto 0) := (others => '0');
    signal merged4_out_roll:    std_logic_vector(NMERGE4-1 downto 0) := (others => '0');
    signal merged4_out_full:    std_logic_vector(NMERGE4-1 downto 0) := (others => '0');

    signal merged_out :        particles(NREGIONS-1 downto 0);
    signal merged_out_valid :  std_logic_vector(NREGIONS-1 downto 0) := (others => '0');
    signal merged_out_roll:    std_logic_vector(NREGIONS-1 downto 0) := (others => '0');

begin

    router : entity work.calo_router 
                port map(ap_clk => ap_clk, 
                             enabled => ap_start,
                             newevent => newevent,
                             links_in => links_in,
                             fifo_in => fifo_in,
                             fifo_in_write => fifo_in_write,
                             fifo_in_roll  => fifo_in_roll);

    gen_fifos: for ireg in NALLFIFOS-1 downto 0 generate
        reg_buffer : entity work.rolling_fifo
                        --generic map(FIFO_INDEX => ireg+1)
                        port map(ap_clk => ap_clk, 
                                 d_in    => fifo_in(ireg),
                                 write_in  => fifo_in_write(ireg),
                                 roll   => fifo_in_roll(ireg),
                                 d_out    => fifo_out(ireg),
                                 valid_out  => fifo_out_valid(ireg),
                                 --dbg_w64 =>  fifo_dbg(ireg),
                                 full  => fifo_out_full(ireg),
                                 roll_out  => fifo_out_roll(ireg)
                             );
        end generate gen_fifos;

    gen_merger2s: for imerge in NMERGE2-1 downto 0 generate
        reg_merger2 : entity work.fifo_merge2_full
                        --generic map(FIFO_INDEX => imerge+1)
                        port map(ap_clk => ap_clk, 
                                 d1_in => fifo_out(imerge*2),
                                 d2_in => fifo_out(imerge*2+1),
                                 d1_valid => fifo_out_valid(imerge*2),
                                 d2_valid => fifo_out_valid(imerge*2+1),
                                 roll     => fifo_out_roll(imerge*2),
                                 full     => merged2_out_full(imerge),
                                 d_out      => merged2_out(imerge),
                                 valid_out  => merged2_out_valid(imerge),
                                 full1      => fifo_out_full(imerge*2),  
                                 full2      => fifo_out_full(imerge*2+1),
                                 --dbg_w64    => merged2_dbg(imerge),
                                 roll_out   => merged2_out_roll(imerge)
                            );
        end generate gen_merger2s;

    gen_merger4s: for imerge in NMERGE4-1 downto 0 generate
        reg_merger4 : entity work.fifo_merge2_full
                        --generic map(FIFO_INDEX => imerge+1)
                        port map(ap_clk => ap_clk, 
                                 d1_in => merged2_out(imerge*2),
                                 d2_in => merged2_out(imerge*2+1),
                                 d1_valid => merged2_out_valid(imerge*2),
                                 d2_valid => merged2_out_valid(imerge*2+1),
                                 roll     => merged2_out_roll(imerge*2),
                                 full     => merged4_out_full(imerge),
                                 d_out      => merged4_out(imerge),
                                 valid_out  => merged4_out_valid(imerge),
                                 full1      => merged2_out_full(imerge*2),  
                                 full2      => merged2_out_full(imerge*2+1),
                                 --dbg_w64    => merged4_dbg(imerge),
                                 roll_out   => merged4_out_roll(imerge)
                            );
        end generate gen_merger4s;

    gen_mergers_s: for isec in NCALOSECTORS-1 downto 0 generate
        delay_0 : process(ap_clk)
            constant iFROM : natural := 5*isec;
            constant iTO : natural := 3*isec;
        begin
            if rising_edge(ap_clk) then
                merged_out(iTO) <= merged4_out(iFROM);
                merged_out_valid(iTO) <= merged4_out_valid(iFROM);
                merged_out_roll(iTO) <= merged4_out_roll(iFROM);
                merged4_out_full(iFROM) <= '0';
            end if;
        end process delay_0;
            
        gen_mergers_12 : for imerge in 0 to 1 generate
            reg_merger : entity work.fifo_merge2
                            port map(ap_clk => ap_clk, 
                                     d1_in => merged4_out(5*isec+2*imerge+1),
                                     d2_in => merged4_out(5*isec+2*imerge+2),
                                     d1_valid => merged4_out_valid(5*isec+2*imerge+1),
                                     d2_valid => merged4_out_valid(5*isec+2*imerge+2),
                                     roll     => merged4_out_roll(3*isec+imerge+1),
                                     --full     => '0',
                                     d_out      => merged_out(3*isec+imerge+1),
                                     valid_out  => merged_out_valid(3*isec+imerge+1),
                                     full1      => merged4_out_full(5*isec+2*imerge+1),  
                                     full2      => merged4_out_full(5*isec+2*imerge+2),
                                     --dbg_w64    => merged_dbg(3*isec+imerge+1),
                                     roll_out   => merged_out_roll(3*isec+imerge+1)
                                );
            end generate gen_mergers_12;
        end generate gen_mergers_s;


    links_in( 0) <= w64_to_particle(calo_in_0_0_V);
    links_in( 1) <= w64_to_particle(calo_in_0_1_V);
    links_in( 2) <= w64_to_particle(calo_in_0_2_V);
    links_in( 3) <= w64_to_particle(calo_in_0_3_V);
    links_in( 4) <= w64_to_particle(calo_in_1_0_V);
    links_in( 5) <= w64_to_particle(calo_in_1_1_V);
    links_in( 6) <= w64_to_particle(calo_in_1_2_V);
    links_in( 7) <= w64_to_particle(calo_in_1_3_V);
    links_in( 8) <= w64_to_particle(calo_in_2_0_V);
    links_in( 9) <= w64_to_particle(calo_in_2_1_V);
    links_in(10) <= w64_to_particle(calo_in_2_2_V);
    links_in(11) <= w64_to_particle(calo_in_2_3_V);

    calo_out_0_V <= particle_to_w64(merged_out(0));
    calo_out_1_V <= particle_to_w64(merged_out(1));
    calo_out_2_V <= particle_to_w64(merged_out(2));
    calo_out_3_V <= particle_to_w64(merged_out(3));
    calo_out_4_V <= particle_to_w64(merged_out(4));
    calo_out_5_V <= particle_to_w64(merged_out(5));
    calo_out_6_V <= particle_to_w64(merged_out(6));
    calo_out_7_V <= particle_to_w64(merged_out(7));
    calo_out_8_V <= particle_to_w64(merged_out(8));
    calo_out_valid_0 <= merged_out_valid(0);
    calo_out_valid_1 <= merged_out_valid(1);
    calo_out_valid_2 <= merged_out_valid(2);
    calo_out_valid_3 <= merged_out_valid(3);
    calo_out_valid_4 <= merged_out_valid(4);
    calo_out_valid_5 <= merged_out_valid(5);
    calo_out_valid_6 <= merged_out_valid(6);
    calo_out_valid_7 <= merged_out_valid(7);
    calo_out_valid_8 <= merged_out_valid(8);

    newevent_out <= merged_out_roll(0);

end Behavioral;
